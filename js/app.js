const windowSize = 25;
const regionLength = 20;
const cutoff = 0.9;

const inputSequenceEl = document.getElementById('input-sequence');
const submitEl = document.getElementById('submit');
const clearEl = document.getElementById('clear');
const resetEl = document.getElementById('reset');

const downloadRegionsEl = document.getElementById('download-regions-csv');
const downloadResiduesEl = document.getElementById('download-residues-csv');

const outputSequenceLengthEl = document.getElementById('output-sequence-length');
const outputSumPWindowsEl = document.getElementById('output-sum-p-windows');
const outputSumPWindows2El = document.getElementById('output-sum-p-windows-2');

const proteinChartEl = document.getElementById('protein-chart');
const proteinChart2El = document.getElementById('protein-pi-q-chart');
const proteinChart3El = document.getElementById('protein-bar');

var protein = null;

const regionsTable = new Tabulator('#regions-table', {
    autoResize: true,
    maxHeight: '300px',
    layout: 'fitColumns',
    columns: [
        {
            title: 'Region number',
            field: 'index',
            sorter: 'number',
        },
        {
            title: 'Region classifier',
            field: 'region',
            sorter: 'string',
        },
        {
            title: 'First residue of region',
            field: 'start',
            sorter: 'number',
        },
        {
            title: 'Last residue of region',
            field: 'end',
            sorter: 'number',
        },
        {
            title: 'Region length',
            field: 'length',
            sorter: 'number',
        },
    ],
});
const residuesTable = new Tabulator('#residues-table', {
    autoResize: true,
    maxHeight: '300px',
    layout: 'fitColumns',
    columns: [
        {
            title: 'Residue number',
            field: 'index',
            sorter: 'number',
        },
        {
            title: 'Amino Acid type',
            field: 'amino_acid',
            sorter: 'string',
        },
        {
            title: 'Residue label',
            field: 'region',
            sorter: 'number',
        },
        {
            title: 'dist_norm',
            field: 'dist_norm',
            sorter: 'number',
        },
        {
            title: 'Residue label π + q',
            field: 'region_pi_q',
            sorter: 'number',
        },
        {
            title: 'dist_norm_pi_q',
            field: 'dist_norm_pi_q',
            sorter: 'number',
        },
    ],
});

const proteinChartLayout = {
    showlegend: true,
    xaxis: {
        title: {
            text: 'Residue Number',
        },
    },
    yaxis: {
        title: {
            text: 'Residue-level dist_norm',
        },
    },
    shapes: [],
};
const proteinChart2Layout = {
    showlegend: true,
    xaxis: {
        title: {
            text: 'Residue Number',
        },
    },
    yaxis: {
        title: {
            text: 'Residue-level dist_norm_π_q',
        },
    },
    shapes: [],
};
const proteinChart3Layout = {
    barmode: 'stack',
    xaxis: {
        title: {
            text: 'Residue Number',
        },
        zeroline: false,
    },
    yaxis: {
        visible: false,
    },
};

function excludeWindow(index, length, windowSize) {
    if (index < windowSize / 2 - 1 ||
        index > length - windowSize / 2) {
        return true;
    }
    return false;
}

function resetSummary() {
    outputSequenceLengthEl.innerText = '';
    outputSumPWindowsEl.innerText = '';
    outputSumPWindows2El.innerText = '';
}
function resetPlots() {
    Plotly.newPlot(proteinChartEl, [], proteinChartLayout);
    Plotly.newPlot(proteinChart2El, [], proteinChart2Layout);
    Plotly.newPlot(proteinChart3El, [], proteinChart3Layout);
}
function resetTables() {
    regionsTable.setData([]);
    residuesTable.setData([]);
}

function displayProteinReuslts() {
    const filteredResidues = protein.residues.filter(x => !excludeWindow(x.index, protein.residues.length, windowSize));
    const regions = protein.regions;

    outputSequenceLengthEl.innerText = protein.sequence.length;
    outputSumPWindowsEl.innerText = filteredResidues
        .map(x => {
            return x.region == 'P' ? x.dist_norm : 0
        })
        .reduce((prev, curr) => prev + curr, 0)
        .toFixed(3);
    outputSumPWindows2El.innerText = filteredResidues
        .map(x => {
            return x.region_pi_q == 'P' ? x.dist_norm_pi_q : 0
        })
        .reduce((prev, curr) => prev + curr, 0)
        .toFixed(3);

    var regionData = [];
    protein.regions.forEach((x, index) => {
        regionData.push({
            index: index + 1,
            region: x.region,
            start: x.start + 1,
            end: x.end + 1,
            length: (x.end - x.start) + 1,
        });
    });
    console.log(regionData);
    regionsTable.setData(regionData);

    var residueData = [];
    protein.residues.forEach((x, index) => {
        residueData.push({
            index: index + 1,
            amino_acid: x.amino_acid,
            region: x.region,
            dist_norm: x.dist_norm.toFixed(3),
            region_pi_q: x.region_pi_q,
            dist_norm_pi_q: x.dist_norm_pi_q.toFixed(3),
        });
    });
    residuesTable.setData(residueData);

    const traces = {
        P: {
            name: 'P labeled residue',
            type: 'scatter',
            mode: 'markers',
            hoverinfo: 'x+y',
            line: {
                color: '#5075B0',
            },
            x: [],
            y: [],
        },
        D: {
            name: 'D labeled residue',
            type: 'scatter',
            mode: 'markers',
            hoverinfo: 'x+y',
            line: {
                color: '#A92217',
            },
            x: [],
            y: [],
        },
        F: {
            name: 'F labeled residue',
            type: 'scatter',
            mode: 'markers',
            hoverinfo: 'x+y',
            line: {
                color: '#000000',
            },
            x: [],
            y: [],
        },
    };

    filteredResidues.forEach((x, index) => {
        traces.P.x.push(index + 1);
        traces.P.y.push(x.region == 'P' ? x.dist_norm : null);
        traces.D.x.push(index + 1);
        traces.D.y.push(x.region == 'D' ? x.dist_norm : null);
        traces.F.x.push(index + 1);
        traces.F.y.push(x.region == 'F' ? x.dist_norm : null);
    });

    const data = [traces.P, traces.D, traces.F];
    Plotly.newPlot(proteinChartEl, data, proteinChartLayout);

    const traces2 = {
        P: {
            name: 'P labeled residue',
            type: 'scatter',
            mode: 'markers',
            hoverinfo: 'x+y',
            line: {
                color: '#5075B0',
            },
            x: [],
            y: [],
        },
        D: {
            name: 'D labeled residue',
            type: 'scatter',
            mode: 'markers',
            hoverinfo: 'x+y',
            line: {
                color: '#A92217',
            },
            x: [],
            y: [],
        },
        F: {
            name: 'F labeled residue',
            type: 'scatter',
            mode: 'markers',
            hoverinfo: 'x+y',
            line: {
                color: '#000000',
            },
            x: [],
            y: [],
        },
    };

    filteredResidues.forEach((x, index) => {
        traces2.P.x.push(index + 1);
        traces2.P.y.push(x.region == 'P' ? x.dist_norm_pi_q : null);
        traces2.D.x.push(index + 1);
        traces2.D.y.push(x.region == 'D' ? x.dist_norm_pi_q : null);
        traces2.F.x.push(index + 1);
        traces2.F.y.push(x.region == 'F' ? x.dist_norm_pi_q : null);
    });

    const data2 = [traces2.P, traces2.D, traces2.F];
    Plotly.newPlot(proteinChart2El, data2, proteinChart2Layout);

    const traceLegendSettings = {
        P: {
            ranked: false,
            rank: 3,
        },
        D: {
            ranked: false,
            rank: 2,
        },
        F: {
            ranked: false,
            rank: 1,
        },
    };

    const traces3 = [];
    for (let i = 0; i < regions.length; i++) {
        if (i == 0 || regions[i].start > (regions[i - 1].end + 1)) {
            var spacerTrace = {
                type: 'bar',
                orientation: 'h',
                width: [0.3],
                x: [],
                y: ['Region'],
                marker: {
                    color: 'transparent',
                },
                showlegend: false,
            };

            if (i == 0) {
                spacerTrace.x.push(regions[i].start);
            } else {
                spacerTrace.x.push(regions[i].start - regions[i - 1].end);
            }

            traces3.push(spacerTrace);
        }


        var trace = {
            type: 'bar',
            orientation: 'h',
            width: [0.3],
            x: [regions[i].end - regions[i].start],
            y: ['Region'],
            marker: {
                color: 'transparent',
            },
            showlegend: false,
        };

        switch (regions[i].region) {
            case 'P': {
                trace.name = 'PS';
                trace.marker.color = '#5075B0';
                if (!traceLegendSettings.P.ranked) {
                    traceLegendSettings.P.ranked = true;
                    trace.showlegend = true;
                    trace.legendrank = traceLegendSettings.P.rank;
                }
            }
                break;
            case 'D': {
                trace.name = 'ID';
                trace.marker.color = '#A92217';
                if (!traceLegendSettings.D.ranked) {
                    traceLegendSettings.D.ranked = true;
                    trace.showlegend = true;
                    trace.legendrank = traceLegendSettings.D.rank;
                }
            }
                break;
            case 'F': {
                trace.name = 'Folded';
                trace.marker.color = '#000000';
                if (!traceLegendSettings.F.ranked) {
                    traceLegendSettings.F.ranked = true;
                    trace.showlegend = true;
                    trace.legendrank = traceLegendSettings.F.rank;
                }
            }
                break;
        }

        traces3.push(trace);
    }

    const data3 = traces3;
    Plotly.newPlot(proteinChart3El, data3, proteinChart3Layout);
}

function processProteinSequence(sequence) {
    var parseResult = parse(sequence, windowSize, regionLength, cutoff);

    if (parseResult.error) {
        console.log(parseResult.message);
    } else {
        protein = {
            sequence: parseResult.data.sequence,
            regions: parseResult.data.regions,
            residues: parseResult.data.residues,
        };
    }
}

function submit() {
    var sequence = inputSequenceEl.value;

    processProteinSequence(sequence);
    displayProteinReuslts();
}
function clear() {
    console.log('clear');
    inputSequenceEl.value = '';
}
function reset() {
    clear();
    resetSummary();
    resetPlots();
    resetTables();

}

reset();

submitEl.addEventListener('click', submit);
clearEl.addEventListener('click', clear);
resetEl.addEventListener('click', reset);
downloadRegionsEl.addEventListener('click', _ => {
    regionsTable.download("csv", "regions.csv");
});
downloadResiduesEl.addEventListener('click', _ => {
    residuesTable.download("csv", "residues.csv");
});
