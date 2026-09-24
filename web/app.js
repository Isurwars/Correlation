/**
 * Correlation WASM — Web Application Logic
 *
 * Handles file upload, coordinates calculation via Web Worker (off-main-thread)
 * with graceful main-thread fallback, and renders interactive distribution function
 * plots on a <canvas> element.
 */

// Global state
let Module = null;
let trajectory = null;
let df = null;
let worker = null;
let useWorker = false;
let cachedResults = null;

// Color palette (Okabe-Ito)
const COLORS = [
    '#E69F00', '#56B4E9', '#009E73', '#F0E442',
    '#0072B2', '#D55E00', '#CC79A7', '#FFFFFF'
];

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
    const dropZone = document.getElementById('drop-zone');
    const fileInput = document.getElementById('file-input');
    const fileInfo = document.getElementById('file-info');
    const controlsPanel = document.getElementById('controls-panel');
    const runBtn = document.getElementById('run-btn');
    const statusEl = document.getElementById('status');
    const resultsPanel = document.getElementById('results-panel');
    const plotSelect = document.getElementById('plot-select');
    const partialSelect = document.getElementById('partial-select');

    // Feature detection for WebAssembly
    if (typeof WebAssembly !== 'object' || typeof WebAssembly.instantiate !== 'function') {
        statusEl.textContent = 'Error: WebAssembly is not supported by your browser.';
        return;
    }

    // Try initializing Web Worker first
    if (window.Worker) {
        try {
            worker = new Worker('worker.js');
            worker.onmessage = handleWorkerMessage;
            worker.onerror = (err) => {
                console.warn('Worker error or blocked by CORS, falling back to main thread:', err);
                useWorker = false;
                initMainThread();
            };
            // Note: worker.js self-inits and will send READY
        } catch (err) {
            console.warn('Could not instantiate worker, falling back to main thread:', err);
            initMainThread();
        }
    } else {
        initMainThread();
    }

    function initMainThread() {
        if (typeof createCorrelationModule === 'function') {
            createCorrelationModule().then(m => {
                Module = m;
                statusEl.textContent = 'Correlation WASM module ready (main thread).';
            }).catch(err => {
                statusEl.textContent = `Error initializing WASM module: ${err.message || err}`;
                console.error('WASM module loading error:', err);
            });
        } else {
            statusEl.textContent = 'WASM module not found — ensure correlation_wasm.js is built.';
        }
    }

    function handleWorkerMessage(e) {
        const { type, message, nFrames, results } = e.data;
        switch (type) {
            case 'READY':
                useWorker = true;
                statusEl.textContent = 'Correlation WASM module ready (Worker thread).';
                break;
            case 'STATUS':
                statusEl.textContent = message;
                break;
            case 'FILE_LOADED':
                controlsPanel.classList.remove('hidden');
                runBtn.disabled = false;
                statusEl.textContent = `Parsed ${nFrames} frame(s). Ready to analyze.`;
                break;
            case 'ANALYSIS_COMPLETE':
                cachedResults = results;
                populateSelectorsFromWorker(results);
                resultsPanel.classList.remove('hidden');
                statusEl.textContent = 'Analysis complete.';
                runBtn.disabled = false;
                renderPlot();
                break;
            case 'ERROR':
                statusEl.textContent = `Error: ${message}`;
                runBtn.disabled = false;
                break;
            default:
                console.warn('Unhandled worker message:', e.data);
        }
    }

    function populateSelectorsFromWorker(results) {
        plotSelect.innerHTML = '';
        const names = Object.keys(results);
        for (const name of names) {
            const opt = document.createElement('option');
            opt.value = name;
            opt.textContent = name;
            plotSelect.appendChild(opt);
        }
    }

    // Drop zone events
    dropZone.addEventListener('click', () => fileInput.click());
    dropZone.addEventListener('dragover', e => {
        e.preventDefault();
        dropZone.classList.add('dragover');
    });
    dropZone.addEventListener('dragleave', () => dropZone.classList.remove('dragover'));
    dropZone.addEventListener('drop', e => {
        e.preventDefault();
        dropZone.classList.remove('dragover');
        if (e.dataTransfer.files.length > 0) handleFile(e.dataTransfer.files[0]);
    });
    fileInput.addEventListener('change', () => {
        if (fileInput.files.length > 0) handleFile(fileInput.files[0]);
    });

    // Run button
    runBtn.addEventListener('click', runAnalysis);

    // Plot selectors
    plotSelect.addEventListener('change', renderPlot);
    partialSelect.addEventListener('change', renderPlot);

    // ---------------------------------------------------------------------------
    // Handle uploaded file
    // ---------------------------------------------------------------------------
    function handleFile(file) {
        fileInfo.classList.remove('hidden');
        fileInfo.textContent = `Loaded: ${file.name} (${(file.size / 1024).toFixed(1)} KB)`;

        const reader = new FileReader();
        reader.onload = () => {
            const data = new Uint8Array(reader.result);
            const strData = new TextDecoder().decode(data);

            if (useWorker && worker) {
                statusEl.textContent = 'Loading and parsing file in worker...';
                worker.postMessage({
                    action: 'LOAD_FILE',
                    payload: { text: strData, filename: file.name }
                });
            } else {
                try {
                    if (!Module) {
                        statusEl.textContent = 'Error: WASM module not ready.';
                        return;
                    }
                    trajectory = Module.readFromBuffer(strData, file.name);
                    const nFrames = trajectory.numFrames();
                    controlsPanel.classList.remove('hidden');
                    runBtn.disabled = false;
                    statusEl.textContent = `Parsed ${nFrames} frame(s). Ready to analyze.`;
                } catch (err) {
                    statusEl.textContent = `Error: ${err.message || err}`;
                }
            }
        };
        reader.readAsArrayBuffer(file);
    }

    // ---------------------------------------------------------------------------
    // Run analysis
    // ---------------------------------------------------------------------------
    function runAnalysis() {
        const rMax = parseFloat(document.getElementById('r-max').value) || 20.0;
        const binWidth = parseFloat(document.getElementById('bin-width').value) || 0.05;
        const calcPad = !!(document.getElementById('calc-pad') && document.getElementById('calc-pad').checked);

        statusEl.textContent = 'Running analysis...';
        runBtn.disabled = true;

        if (useWorker && worker) {
            worker.postMessage({
                action: 'RUN_ANALYSIS',
                payload: { rMax, binWidth, calcPad }
            });
        } else {
            if (!trajectory || !Module) return;
            setTimeout(() => {
                try {
                    df = new Module.DistributionFunctions(trajectory, 0.0, []);
                    df.calculateRDF(rMax, binWidth);
                    if (calcPad) {
                        df.calculatePAD(0.5);
                    }

                    const histNames = df.getAvailableHistograms();
                    plotSelect.innerHTML = '';
                    for (let i = 0; i < histNames.size(); i++) {
                        const opt = document.createElement('option');
                        opt.value = histNames.get(i);
                        opt.textContent = histNames.get(i);
                        plotSelect.appendChild(opt);
                    }

                    resultsPanel.classList.remove('hidden');
                    statusEl.textContent = 'Analysis complete.';
                    renderPlot();
                } catch (err) {
                    statusEl.textContent = `Error: ${err.message || err}`;
                }
                runBtn.disabled = false;
            }, 50);
        }
    }

    // ---------------------------------------------------------------------------
    // Render plot on <canvas>
    // ---------------------------------------------------------------------------
    function renderPlot() {
        const histName = plotSelect.value;
        if (!histName) return;

        let bins = null;
        let keys = [];
        let title = histName;
        let xLabel = '';
        let yLabel = '';
        let getPartialData = null;

        if (useWorker && cachedResults) {
            const histData = cachedResults[histName];
            if (!histData) return;
            bins = histData.bins;
            keys = histData.keys;
            title = histData.title || histName;
            xLabel = histData.xLabel || '';
            yLabel = histData.yLabel || '';
            getPartialData = (k) => histData.partials[k];
        } else if (df) {
            const hist = df.getHistogram(histName);
            bins = Array.from(hist.getBins());
            keys = Array.from(hist.getPartialKeys());
            title = hist.title || histName;
            xLabel = hist.xLabel || '';
            yLabel = hist.yLabel || '';
            getPartialData = (k) => {
                const p = hist.getPartial(k);
                return p ? Array.from(p) : null;
            };
        } else {
            return;
        }

        // Populate partial selector
        const currentPartial = partialSelect.value;
        partialSelect.innerHTML = '';
        for (let i = 0; i < keys.length; i++) {
            const opt = document.createElement('option');
            opt.value = keys[i];
            opt.textContent = keys[i];
            if (keys[i] === currentPartial) opt.selected = true;
            partialSelect.appendChild(opt);
        }
        if (!partialSelect.value && keys.length > 0) {
            for (let i = 0; i < keys.length; i++) {
                if (keys[i] === 'Total') { partialSelect.value = 'Total'; break; }
            }
            if (!partialSelect.value) partialSelect.value = keys[0];
        }

        const partialKey = partialSelect.value;
        const ys = getPartialData(partialKey);
        if (!ys) return;

        drawChart(bins, ys, title, xLabel, yLabel, partialKey);
    }

    // ---------------------------------------------------------------------------
    // Canvas chart renderer
    // ---------------------------------------------------------------------------
    function drawChart(xs, ys, title, xLabel, yLabel, seriesLabel) {
        const canvas = document.getElementById('plot-canvas');
        const ctx = canvas.getContext('2d');
        const W = canvas.width;
        const H = canvas.height;
        const pad = { top: 50, right: 30, bottom: 60, left: 80 };

        ctx.clearRect(0, 0, W, H);
        ctx.fillStyle = '#1a1a2e';
        ctx.fillRect(0, 0, W, H);

        const n = Math.min(xs.length, ys.length);
        if (n === 0) return;

        let xMin = xs[0], xMax = xs[n - 1];
        let yMin = 0, yMax = 0;
        for (let i = 0; i < n; i++) {
            if (ys[i] > yMax) yMax = ys[i];
            if (ys[i] < yMin) yMin = ys[i];
        }
        yMax *= 1.05;

        const px = v => pad.left + (v - xMin) / (xMax - xMin) * (W - pad.left - pad.right);
        const py = v => pad.top + (1 - (v - yMin) / (yMax - yMin)) * (H - pad.top - pad.bottom);

        // Grid
        ctx.strokeStyle = 'rgba(255,255,255,0.06)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 5; i++) {
            const y = pad.top + i / 5 * (H - pad.top - pad.bottom);
            ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(W - pad.right, y); ctx.stroke();
        }

        // Axes
        ctx.strokeStyle = '#cdd6f4';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(pad.left, pad.top);
        ctx.lineTo(pad.left, H - pad.bottom);
        ctx.lineTo(W - pad.right, H - pad.bottom);
        ctx.stroke();

        // Data line
        ctx.strokeStyle = COLORS[0];
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < n; i++) {
            const x = px(xs[i]);
            const y = py(ys[i]);
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        }
        ctx.stroke();

        // Labels
        ctx.fillStyle = '#a6adc8';
        ctx.font = '14px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(xLabel || '', (pad.left + W - pad.right) / 2, H - 15);

        ctx.save();
        ctx.translate(20, (pad.top + H - pad.bottom) / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText(yLabel || '', 0, 0);
        ctx.restore();

        // Title
        ctx.font = '16px Inter, sans-serif';
        ctx.fillStyle = '#e0e0e0';
        ctx.fillText(title, W / 2, 30);

        // Legend
        ctx.font = '12px Inter, sans-serif';
        ctx.fillStyle = COLORS[0];
        ctx.textAlign = 'left';
        ctx.fillRect(W - pad.right - 100, pad.top + 5, 15, 3);
        ctx.fillStyle = '#a6adc8';
        ctx.fillText(seriesLabel, W - pad.right - 80, pad.top + 12);
    }
});
