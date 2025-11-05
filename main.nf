println "input_xlsx: ${params.input_xlsx}"
println "input_pretty_names: ${params.input_pretty_names}"

workflow {
  Channel
    .fromPath(params.input_xlsx)
    .set { input_xlsx_ch }

  Channel
    .fromPath(params.input_pretty_names)
    .set { input_pretty_names_ch }

  load_data(input_xlsx_ch, input_pretty_names_ch)
}


process load_data {
  input:
  path input_xlsx_ch
  path input_pretty_names_ch

  output:
  path 'processed_data/data_df_pre_scaling.rds'
  path 'processed_data/data_df.rds'
  path 'processed_data/data_df_renamed.rds'
  path 'processed_data/variables.rds'
  path 'processed_data/data_df_complete.rds'

  script:
  """

  Rscript -e 'source("${projectDir}/00-universal-dependencies.R"); source("${projectDir}/01-load-data.R")'
  mkdir -p processed_data
  cp ${projectDir}/${params.outdir}/processed_data/data_df_pre_scaling.rds processed_data/data_df_pre_scaling.rds
  cp ${projectDir}/${params.outdir}/processed_data/data_df.rds processed_data/data_df.rds
  cp ${projectDir}/${params.outdir}/processed_data/data_df_renamed.rds processed_data/data_df_renamed.rds
  cp ${projectDir}/${params.outdir}/processed_data/variables.rds processed_data/variables.rds
  cp ${projectDir}/${params.outdir}/processed_data/data_df_complete.rds processed_data/data_df_complete.rds

  """
}

