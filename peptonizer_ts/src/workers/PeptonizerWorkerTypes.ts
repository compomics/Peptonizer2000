// Define a specific type of inputs that are expected for each task that can be performed by this worker.
enum WorkerTask {
    PERFORM_TAXA_WEIGHING,
    GENERATE_GRAPH,
    EXECUTE_PEPGM,
    CLUSTER_TAXA,
    FIND_BEST_PARAMETERS
}

interface PerformTaxaWeighingTaskData {
    peptidesScores: Map<string, number>;
    peptidesCounts: Map<string, number>;
}

interface GenerateGraphTaskData {
    taxaWeightsCsv: string;
}

interface ExecutePepgmTaskData {
    graphXml: string,
    alpha: number,
    beta: number,
    prior: number
}

interface ClusterTaxaTaskData {

}

interface FindBestParametersTaskData {

}

type SpecificInputEventData =
    { task: WorkerTask.PERFORM_TAXA_WEIGHING, input: PerformTaxaWeighingTaskData } |
    { task: WorkerTask.GENERATE_GRAPH, input: GenerateGraphTaskData } |
    { task: WorkerTask.EXECUTE_PEPGM, input: ExecutePepgmTaskData } |
    { task: WorkerTask.CLUSTER_TAXA, input: ClusterTaxaTaskData } |
    { task: WorkerTask.FIND_BEST_PARAMETERS, input: FindBestParametersTaskData };

type CommonInputEventData = { workerId: number };

type InputEventData = SpecificInputEventData & CommonInputEventData;

interface PerformTaxaWeighingTaskResult {
    taxaWeightsCsv: string
}

interface GenerateGraphTaskDataResult {
    graphXml: string
}

interface ExecutePepgmTaskDataResult {
    taxonScoresJson: string
}

interface PepgmProgressUpdate {
    progressType: "graph" | "max_residual" | "iteration",
    currentValue: number,
    maxValue: number
}

interface ClusterTaxaTaskDataResult {

}

interface FindBestParametersTaskDataResult {

}

enum ResultType {
    SUCCESSFUL,
    PROGRESS,
    FAILED
}

type SpecificOutputEventData = { resultType: ResultType.SUCCESSFUL } & (
    { task: WorkerTask.PERFORM_TAXA_WEIGHING, output: PerformTaxaWeighingTaskResult } |
    { task: WorkerTask.GENERATE_GRAPH, output: GenerateGraphTaskDataResult } |
    { task: WorkerTask.EXECUTE_PEPGM, output: ExecutePepgmTaskDataResult } |
    { task: WorkerTask.CLUSTER_TAXA, output: ClusterTaxaTaskDataResult } |
    { task: WorkerTask.FIND_BEST_PARAMETERS, output: FindBestParametersTaskDataResult });

type CommonOutputEventData = { workerId: number };

type ErrorOutputEvent = { resultType: ResultType.FAILED, error: string };
type ProgressOutputEvent = { resultType: ResultType.PROGRESS, task: WorkerTask.EXECUTE_PEPGM, progressUpdate: PepgmProgressUpdate };

type OutputEventData = (SpecificOutputEventData | ErrorOutputEvent | ProgressOutputEvent) & CommonOutputEventData;

export {
    WorkerTask,
    ResultType
};

export type {
    PerformTaxaWeighingTaskData,
    GenerateGraphTaskData,
    ExecutePepgmTaskData,
    ClusterTaxaTaskData,
    FindBestParametersTaskData,
    SpecificInputEventData,
    InputEventData,
    PerformTaxaWeighingTaskResult,
    GenerateGraphTaskDataResult,
    ExecutePepgmTaskDataResult,
    ClusterTaxaTaskDataResult,
    FindBestParametersTaskDataResult,
    OutputEventData,
    PepgmProgressUpdate,
    ProgressOutputEvent
};
