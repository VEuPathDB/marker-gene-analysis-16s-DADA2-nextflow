#!/usr/bin/env bash

mkdir filtered
Rscript /usr/bin/filterFastqs.R \
  --fastqsInDir . \
  --fastqsOutDir ./filtered \
  --isPaired $params.isPaired \
  --trimLeft $params.trimLeft \
  --trimLeftR $params.trimLeftR \
  --truncLen $params.truncLen \
  --truncLenR $params.truncLenR \
  --maxLen $params.maxLen \
  --platform $params.platform
