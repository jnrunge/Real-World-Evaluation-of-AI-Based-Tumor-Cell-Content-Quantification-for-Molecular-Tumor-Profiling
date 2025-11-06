params.input_xlsx = "input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx"
params.input_pretty_names = "input/variable_pretty_names.csv"
params.input_merged_data = "input/raw/merged_data.csv"

workflow {
  input_xlsx_ch = Channel.fromPath(params.input_xlsx)
  input_pretty_names_ch = Channel.fromPath(params.input_pretty_names)
  input_merged_data_ch = Channel.fromPath(params.input_merged_data)
  rscript1_ch = Channel.fromPath("00-universal-dependencies.R")
  rscript2_ch = Channel.fromPath("01-load-data.R")
  rscript3_ch = Channel.fromPath("02-univariate-analyses.R")
  rscript4_ch = Channel.fromPath("02b-tobi-inspired-plot.R")
  rscript5_ch = Channel.fromPath("02c-other-univ-plots.R")
  rscript6_ch = Channel.fromPath("03-bootstrapped-model-inclusion.R")
  rscript7_ch = Channel.fromPath("04-best-models.R")
  rscript7a_ch = Channel.fromPath("04a-plot-models.R")
  rscript8_ch = Channel.fromPath("05-clustering.R")
  function1_ch = Channel.fromPath("functions/manual-stepwise.R")
  function2_ch = Channel.fromPath("functions/bootstrapping_models.R")
  function3_ch = Channel.fromPath("functions/best_models.R")

  data_outputs = load_data(input_xlsx_ch, input_pretty_names_ch, input_merged_data_ch, rscript1_ch, rscript2_ch)
  univariate_outputs = univariate_analysis(input_xlsx_ch, input_pretty_names_ch, data_outputs, rscript1_ch, rscript2_ch, rscript3_ch, rscript4_ch, rscript5_ch)
  bootstrapped_outputs = bootstrapped_models(input_xlsx_ch, input_pretty_names_ch, data_outputs, rscript1_ch, rscript2_ch, rscript6_ch, function1_ch, function2_ch)
  best_models_output = best_models(input_pretty_names_ch, data_outputs, rscript1_ch, rscript2_ch, rscript7_ch, rscript7a_ch, function3_ch, function1_ch, input_xlsx_ch)
  clustering_outputs = clustering(input_xlsx_ch, input_pretty_names_ch, input_merged_data_ch, best_models_output, rscript1_ch, rscript2_ch, rscript8_ch)
}

process load_data {
    publishDir "${projectDir}", mode: 'copy', pattern: 'output/**'

  input:
  path input_xlsx_ch, stageAs: 'input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx'
  path input_pretty_names_ch, stageAs: 'input/variable_pretty_names.csv'
  path input_merged_data_ch, stageAs: 'input/raw/merged_data.csv'
  path rscript1_ch
  path rscript2_ch

  output:
  path 'output/processed_data/data_df_pre_scaling.rds'
  path 'output/processed_data/data_df.rds'
  path 'output/processed_data/data_df_renamed.rds'
  path 'output/processed_data/variables.rds'
  path 'output/processed_data/data_df_complete.rds'
  path 'output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds'
  path 'output/outliers/*'
  path 'output/distributions/*'


  script:
  """
  Rscript -e 'source("${rscript1_ch}"); source("${rscript2_ch}")'
  """
}

process univariate_analysis {
    publishDir "${projectDir}", mode: 'copy', pattern: 'output/**'

  input:
  path input_xlsx_ch, stageAs: 'input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx'
  path input_pretty_names_ch, stageAs: 'input/variable_pretty_names.csv'
  path 'output/processed_data/data_df_pre_scaling.rds'
  path 'output/processed_data/data_df.rds'
  path 'output/processed_data/data_df_renamed.rds'
  path 'output/processed_data/variables.rds'
  path 'output/processed_data/data_df_complete.rds'
  path 'output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds'
  path 'output/outliers/*'
  path 'output/distributions/*'
  path rscript1_ch
  path rscript2_ch
  path rscript3_ch
  path rscript4_ch
  path rscript5_ch

  output:
  path 'output/univar/*'
  path 'output/tobi_plots/*'

  script:
  """
  Rscript -e 'source("${rscript1_ch}"); source("${rscript2_ch}"); source("${rscript3_ch}"); source("${rscript4_ch}"); source("${rscript5_ch}")'
  """
}

process bootstrapped_models{
  publishDir "${projectDir}", mode: 'copy', pattern: 'output/**'
maxForks 1
    input:
  path input_xlsx_ch, stageAs: 'input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx'
  path input_pretty_names_ch, stageAs: 'input/variable_pretty_names.csv'
  path 'output/processed_data/data_df_pre_scaling.rds'
  path 'output/processed_data/data_df.rds'
  path 'output/processed_data/data_df_renamed.rds'
  path 'output/processed_data/variables.rds'
  path 'output/processed_data/data_df_complete.rds'
  path 'output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds'
  path 'output/outliers/*'
  path 'output/distributions/*'
  path rscript1_ch
  path rscript2_ch
  path rscript6_ch
  path function1_ch, stageAs: 'functions/manual-stepwise.R'
  path function2_ch, stageAs: 'functions/bootstrapping_models.R'

  output:
  path 'output/processed_data/bootstrapped_models_forward_*.rds'
  path 'output/bootstrapped_models/*'

  script:
  """
  Rscript -e 'source("${rscript1_ch}"); source("${rscript2_ch}"); source("${rscript6_ch}")'
  """
}


process best_models{
    publishDir "${projectDir}", mode: 'copy', pattern: 'output/**'

  maxForks 1
  input:
  path input_pretty_names_ch, stageAs: 'input/variable_pretty_names.csv'
  path 'output/processed_data/data_df_pre_scaling.rds'
  path 'output/processed_data/data_df.rds'
  path 'output/processed_data/data_df_renamed.rds'
  path 'output/processed_data/variables.rds'
  path 'output/processed_data/data_df_complete.rds'
  path 'output/processed_data/data_df_pre_scaling_NAd_sampletypes.rds'
  path 'output/outliers/*'
  path 'output/distributions/*'
  path rscript1_ch
  path rscript2_ch
  path rscript7_ch
  path rscript7a_ch
  path function3_ch, stageAs: 'functions/best_models.R'
  path function1_ch, stageAs: 'functions/manual-stepwise.R'
  path input_xlsx_ch, stageAs: 'input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx'



  output:
  path 'output/processed_data/best_models_classic_no_forced_interactions.rds'
  path 'output/model_plots/*'

  script:
  """
  Rscript -e 'source("${rscript1_ch}"); source("${rscript2_ch}"); source("${rscript7_ch}"); source("${rscript7a_ch}")'
  """
}

process clustering{
    publishDir "${projectDir}", mode: 'copy', pattern: 'output/**'

  maxForks 1
  
  input:
  path input_xlsx_ch, stageAs: 'input/PathAI_Dataset_300_cases_mg_27.10.2025.xlsx'
  path input_pretty_names_ch, stageAs: 'input/variable_pretty_names.csv'
  path input_merged_data_ch, stageAs: 'input/raw/merged_data.csv'
  path 'output/processed_data/best_models_classic_no_forced_interactions.rds'
  path 'output/model_plots/*'
  path rscript1_ch
  path rscript2_ch
  path rscript8_ch

  output:
  path 'output/clustering/*'

  script:
  """
  Rscript -e 'source("${rscript1_ch}"); source("${rscript2_ch}"); source("${rscript8_ch}")'
  """
}

