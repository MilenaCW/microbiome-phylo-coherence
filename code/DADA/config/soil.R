# DADA2 config for soil dataset.

list(
  dataset = "soil",
  input = list(
    directory = "./soil/data/downloaded_data/16S"
  ),
  output = list(
    directory = "./soil/data/processed_data/16S/DADA"
  ),
  patterns = list(
    forward = "_R1.fastq.gz",
    reverse = "_R2.fastq.gz"
  ),
  filtering = list(
    truncQ = 2,
    trimLeft = 20,
    truncLen = c(250, 240),
    maxN = 0,
    maxEE = c(3, 3)
  ),
  sequence_table = list(
    length_range = 375:475,
    expected_length = 425
  ),
  error_model = list(
    method = "standard"
  )
)
