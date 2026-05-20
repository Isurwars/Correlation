/**
 * Correlation WASM — Web Application Logic
 *
 * Handles file upload, invokes the Emscripten-compiled Correlation module,
 * and renders distribution function plots on a <canvas> element.
 */

// Global state
let Module = null;
let trajectory = null;
let df = null;

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

    // Wait for Emscripten module.
    if (typeof createCorrelationModule === 'function') {
        createCorrelationModule().then(m => {
            Module = m;
            statusEl.textContent = 'WASM module loaded.';
        });
    } else {
        statusEl.textContent = 'WASM module not found — ensure correlation_wasm.js is built.';
    }

    // Drop zone events.
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

    // Run button.
    runBtn.addEventListener('click', runAnalysis);

    // Plot selectors.
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
            try {
                if (!Module) {
                    statusEl.textContent = 'Error: WASM module not ready.';
                    return;
                }
                const data = new Uint8Array(reader.result);
                const strData = new TextDecoder().decode(data);
                trajectory = Module.readFromBuffer(strData, file.name);
                const nFrames = trajectory.numFrames();

                controlsPanel.classList.remove('hidden');
                runBtn.disabled = false;
                statusEl.textContent = `Parsed ${nFrames} frame(s). Ready to analyze.`;
            } catch (err) {
                statusEl.textContent = `Error: ${err.message || err}`;
            }
        };
        reader.readAsArrayBuffer(file);
    }

    // ---------------------------------------------------------------------------
    // Run analysis
    // ---------------------------------------------------------------------------
    function runAnalysis() {
        if (!trajectory || !Module) return;
        const rMax = parseFloat(document.getElementById('r-max').value) || 20.0;
        const binWidth = parseFloat(document.getElementById('bin-width').value) || 0.05;

        statusEl.textContent = 'Running analysis...';
        runBtn.disabled = true;

        setTimeout(() => {
            try {
                // Use last frame for analysis.
                const frames = trajectory.numFrames();
                // We construct DF on the last frame.
                // Note: API matches the embind registration.
                df = new Module.DistributionFunctions(
                    trajectory, // will internally get last frame
                    0.0,
                    []
                );
                df.calculateRDF(rMax, binWidth);

                // Populate plot selector.
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

    // ---------------------------------------------------------------------------
    // Render plot on <canvas>
    // ---------------------------------------------------------------------------
    function renderPlot() {
        if (!df) return;
        const histName = plotSelect.value;
        if (!histName) return;

        const hist = df.getHistogram(histName);
        const bins = hist.getBins();
        const keys = hist.getPartialKeys();

        // Populate partial selector.
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
            // Prefer "Total" if available.
            for (let i = 0; i < keys.length; i++) {
                if (keys[i] === 'Total') { partialSelect.value = 'Total'; break; }
            }
            if (!partialSelect.value) partialSelect.value = keys[0];
        }

        const partialKey = partialSelect.value;
        const ys = hist.getPartial(partialKey);
        if (!ys) return;

        drawChart(bins, ys, hist.title || histName, hist.xLabel, hist.yLabel, partialKey);
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

        // Grid.
        ctx.strokeStyle = 'rgba(255,255,255,0.06)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 5; i++) {
            const y = pad.top + i / 5 * (H - pad.top - pad.bottom);
            ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(W - pad.right, y); ctx.stroke();
        }

        // Axes.
        ctx.strokeStyle = '#cdd6f4';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(pad.left, pad.top);
        ctx.lineTo(pad.left, H - pad.bottom);
        ctx.lineTo(W - pad.right, H - pad.bottom);
        ctx.stroke();

        // Data line.
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

        // Labels.
        ctx.fillStyle = '#a6adc8';
        ctx.font = '14px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(xLabel || '', (pad.left + W - pad.right) / 2, H - 15);

        ctx.save();
        ctx.translate(20, (pad.top + H - pad.bottom) / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText(yLabel || '', 0, 0);
        ctx.restore();

        // Title.
        ctx.font = '16px Inter, sans-serif';
        ctx.fillStyle = '#e0e0e0';
        ctx.fillText(title, W / 2, 30);

        // Legend.
        ctx.font = '12px Inter, sans-serif';
        ctx.fillStyle = COLORS[0];
        ctx.textAlign = 'left';
        ctx.fillRect(W - pad.right - 100, pad.top + 5, 15, 3);
        ctx.fillStyle = '#a6adc8';
        ctx.fillText(seriesLabel, W - pad.right - 80, pad.top + 12);
    }
});
