import {loadPyodide, PyodideInterface} from 'pyodide';
import peptonizerWhlBase64 from "./lib/peptonizer-0.1-py3-none-any.base64.whl?raw";
import {
    ClusterTaxaTaskDataResult,
    ExecutePepgmTaskData,
    ExecutePepgmTaskDataResult,
    FindBestParametersTaskDataResult,
    GenerateGraphTaskData,
    GenerateGraphTaskDataResult,
    InputEventData,
    OutputEventData,
    PerformTaxaWeighingTaskData,
    PerformTaxaWeighingTaskResult,
    ResultType,
    WorkerTask
} from "./PeptonizerWorkerTypes.ts";

import performTaxaWeighingPythonCode from "./lib/perform_taxa_weighing.py?raw";
import generateGraphPythonCode from "./lib/generate_pepgm_graph.py?raw";
import executePepgmPythonCode from "./lib/execute_pepgm.py?raw";

interface DedicatedWorkerGlobalScope {
    pyodide: PyodideInterface;
    postMessage: (message: OutputEventData) => void;
    submitPepgmProgress: (progressType: "graph" | "max_residual" | "iteration", currentValue: number, maxValue: number, workerId: number) => void;
}

declare const self: DedicatedWorkerGlobalScope & typeof globalThis;

async function loadPyodideAndPackages(): Promise<void> {
    self.pyodide = await loadPyodide({
        indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.3/full/'
    });
    // Load all packages into the Pyodide runtime environment that are required by the Peptonizer
    await self.pyodide.loadPackage([
        'numpy',
        'scipy',
        'networkx',
        'pandas',
        'micropip',
        'requests',
        'openssl'
    ]);
    // Use the imported .whl file URL directly with micropip
    await self.pyodide.runPythonAsync(`
        import base64
        from pathlib import Path
        
        import micropip
        
        await micropip.install('rbo')

        # Decode base64 string to binary and write to a temporary file
        wheel_data = "${peptonizerWhlBase64}"
        wheel_binary = base64.b64decode(wheel_data)
        
        # Define a temporary path for the .whl file
        wheel_path = Path("/tmp/peptonizer-0.1-py3-none-any.whl")
        wheel_path.write_bytes(wheel_binary)

        # Install the wheel package
        await micropip.install("emfs:///tmp/peptonizer-0.1-py3-none-any.whl")

        # Clean up by deleting the temporary file
        wheel_path.unlink()
    `);
}

async function performTaxaWeighing(data: PerformTaxaWeighingTaskData): Promise<PerformTaxaWeighingTaskResult> {
    // Set inputs for the Python code
    self.pyodide.globals.set('peptides_scores', data.peptidesScores);
    self.pyodide.globals.set('peptides_counts', data.peptidesCounts);

    // Fetch the Python code and execute it with Pyodide
    const csvString = await self.pyodide.runPythonAsync(performTaxaWeighingPythonCode);

    console.log("Finished taxa weighing")

    return {
        taxaWeightsCsv: csvString
    };
}

async function generateGraph(data: GenerateGraphTaskData): Promise<GenerateGraphTaskDataResult> {
    self.pyodide.globals.set('taxa_weights_csv', data.taxaWeightsCsv);

    const graphXml = await self.pyodide.runPythonAsync(generateGraphPythonCode);

    return {
        graphXml
    };
}

async function executePepgm(data: ExecutePepgmTaskData, workerId: number): Promise<ExecutePepgmTaskDataResult> {
    self.pyodide.globals.set('graph', data.graphXml);
    self.pyodide.globals.set('alpha', data.alpha);
    self.pyodide.globals.set('beta', data.beta);
    self.pyodide.globals.set('prior', data.prior);
    self.pyodide.globals.set('worker_id', workerId);

    const taxonScoresJson = await self.pyodide.runPythonAsync(executePepgmPythonCode);

    console.log(taxonScoresJson);

    return {
        taxonScoresJson
    };
}

// async function clusterTaxa() {
//
// }
//
// async function findBestParameters() {
//
// }

self.submitPepgmProgress = function(
    progressType: "graph" | "max_residual" | "iteration",
    currentValue: number,
    maxValue: number,
    workerId: number
) {
    const resultMessage: OutputEventData = {
        resultType: ResultType.PROGRESS,
        task: WorkerTask.EXECUTE_PEPGM,
        workerId: workerId,
        progressUpdate: {
            progressType,
            currentValue,
            maxValue
        }
    }

    self.postMessage(resultMessage);
}

let pyodideReadyPromise: Promise<void> = loadPyodideAndPackages();

self.onmessage = async (event: MessageEvent<InputEventData>): Promise<void> => {
    try {
        // Make sure loading is done
        await pyodideReadyPromise;

        // Destructure the data from the event
        const eventData = event.data;

        let output: PerformTaxaWeighingTaskResult | GenerateGraphTaskDataResult | ExecutePepgmTaskDataResult | ClusterTaxaTaskDataResult | FindBestParametersTaskDataResult | undefined;

        if (eventData.task === WorkerTask.PERFORM_TAXA_WEIGHING) {
            output = await performTaxaWeighing(eventData.input);
        } else if (eventData.task === WorkerTask.GENERATE_GRAPH) {
            output = await generateGraph(eventData.input);
        } else if (eventData.task === WorkerTask.EXECUTE_PEPGM) {
            output = await executePepgm(eventData.input, eventData.workerId);
        } else if (eventData.task === WorkerTask.CLUSTER_TAXA) {

        } else if (eventData.task === WorkerTask.FIND_BEST_PARAMETERS) {

        } else {
            throw new Error("Unknown task type passed to worker!");
        }

        if (!output) {
            throw new Error("No valid output was generated by worker!");
        }

        self.postMessage({
            resultType: ResultType.SUCCESSFUL,
            workerId: eventData.workerId,
            task: eventData.task,
            output
        });
    } catch (error: any) {
        console.error(error);

        self.postMessage({
            resultType: ResultType.FAILED,
            workerId: event.data.workerId,
            error: error.toString()
        });
    }
};
