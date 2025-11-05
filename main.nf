process load_data {
  input:
  path input_xlsx from params.input_xlsx
  path input_pretty_names from params.input_pretty_names

  output:
  path 'output/processed_data/data_df_pre_scaling.rds'
  path "output/processed_data/data_df.rds"
  path "output/processed_data/data_df_renamed.rds"
  path "output/processed_data/variables.rds"
  path "output/processed_data/data_df_complete.rds"

  script:
  """
  Rscript -e 'source("00-universal-dependencies.R"); source("01-load-data.R")'
  """
}

workflow {
  load_data()
}