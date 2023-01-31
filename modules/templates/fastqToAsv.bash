#!/usr/bin/env bash

set -euo pipefail
Rscript /usr/bin/fastqToAsv.R  \
  --fastqsInDir .  \
  --errorsRdsPath ./$errorFile \
  --outRdsPath ./featureTable.rds \
  --isPaired $params.isPaired \
  --platform $params.platform \
  --mergeTechReps $params.mergeTechReps
