clean:
	rm -rf work/
	rm -rf .nextflow.log*
	cd tests && rm -rf work && rm -rf .nextflow.log* && rm -rf .nextflow/

lint:
	nf-core lint .

test: clean lint
	cd tests && ./run_test.sh -d hadrieng/mag:0.1.0dev -t ../data

docker:
	docker rmi --force hadrieng/mag:0.1.0dev
	docker build -t hadrieng/mag:0.1.0dev .
