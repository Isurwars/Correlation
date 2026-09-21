/**
 * Correlation WASM — Dedicated Web Worker
 *
 * Runs file parsing and distribution function calculations in an isolated background thread
 * to guarantee that the browser UI never freezes or stutters during heavy computation.
 */

let Module = null;
let trajectory = null;
let df = null;

// Load Emscripten WASM glue code
try {
    importScripts('correlation_wasm.js');
} catch (err) {
    self.postMessage({ type: 'ERROR', message: 'Failed to import correlation_wasm.js: ' + err.message });
}

function initWasm() {
    if (typeof createCorrelationModule === 'function') {
        createCorrelationModule().then(m => {
            Module = m;
            self.postMessage({ type: 'READY' });
        }).catch(err => {
            self.postMessage({ type: 'ERROR', message: 'Failed to initialize WASM module: ' + (err.message || err) });
        });
    } else {
        self.postMessage({ type: 'ERROR', message: 'createCorrelationModule is not defined.' });
    }
}

self.onmessage = function(e) {
    const { action, payload } = e.data;

    switch (action) {
        case 'INIT':
            initWasm();
            break;

        case 'LOAD_FILE':
            if (!Module) {
                self.postMessage({ type: 'ERROR', message: 'WASM module not ready.' });
                return;
            }
            try {
                const { text, filename } = payload;
                trajectory = Module.readFromBuffer(text, filename);
                const nFrames = trajectory.numFrames();
                self.postMessage({ type: 'FILE_LOADED', nFrames });
            } catch (err) {
                self.postMessage({ type: 'ERROR', message: 'Failed to load file: ' + (err.message || err) });
            }
            break;

        case 'RUN_ANALYSIS':
            if (!Module || !trajectory) {
                self.postMessage({ type: 'ERROR', message: 'No trajectory loaded.' });
                return;
            }
            try {
                const { rMax, binWidth, calcPad } = payload;
                self.postMessage({ type: 'STATUS', message: 'Initializing distribution functions...' });

                df = new Module.DistributionFunctions(trajectory, 0.0, []);
                self.postMessage({ type: 'STATUS', message: 'Calculating Radial Distribution Function (RDF)...' });
                df.calculateRDF(rMax, binWidth);

                if (calcPad) {
                    self.postMessage({ type: 'STATUS', message: 'Calculating Plane Angle Distribution (PAD)...' });
                    df.calculatePAD(0.5);
                }

                // Extract all histogram data for transferable posting
                const histNames = df.getAvailableHistograms();
                const results = {};

                for (let i = 0; i < histNames.size(); i++) {
                    const name = histNames.get(i);
                    const hist = df.getHistogram(name);
                    const bins = Array.from(hist.getBins());
                    const keys = hist.getPartialKeys();
                    const partials = {};

                    for (let k = 0; k < keys.length; k++) {
                        const key = keys[k];
                        const partData = hist.getPartial(key);
                        if (partData) {
                            partials[key] = Array.from(partData);
                        }
                    }

                    results[name] = {
                        title: hist.title,
                        xLabel: hist.xLabel,
                        yLabel: hist.yLabel,
                        bins,
                        partials,
                        keys: Array.from(keys)
                    };
                }

                self.postMessage({ type: 'ANALYSIS_COMPLETE', results });
            } catch (err) {
                self.postMessage({ type: 'ERROR', message: 'Analysis failed: ' + (err.message || err) });
            }
            break;

        default:
            self.postMessage({ type: 'ERROR', message: 'Unknown action: ' + action });
    }
};

// Auto-trigger initialization on worker startup
initWasm();
