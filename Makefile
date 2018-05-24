clean:
	rm -rf work/
	rm -rf .nextflow.log*
	cd tests && rm -rf work && rm -rf .nextflow.log* && rm -rf .nextflow/
