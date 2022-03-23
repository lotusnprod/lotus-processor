.PHONY: gathering-validation
.PRECIOUS: %.zip %.gz %.tsv %.csv

gathering-validation: ${VALIDATION_PATH}/manuallyValidated.tsv.gz
${VALIDATION_PATH}/manuallyValidated.tsv.gz:
	@[ ! -f ${VALIDATION_PATH} ] && mkdir -p ${VALIDATION_PATH} && echo "Downloading" && wget -nv "https://osf.io/${OSF_VALIDATION}/download" -O ${PWD}/data/validation.zip && unzip -o ${PWD}/data/validation.zip -d ${PWD}/data/ && rm ${PWD}/data/validation.zip && echo "Done" || true
