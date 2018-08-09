clean:
	rm -rf work/
	rm -rf .nextflow.log*
	cd tests && rm -rf work && rm -rf .nextflow.log* && rm -rf .nextflow/

test:
	cd tests && ./run_test.sh -d hadrieng/mag:0.1.0 -t ../data

docker:
	docker rmi --force hadrieng/mag:0.1.0
	docker build -t hadrieng/mag:0.1.0 .
