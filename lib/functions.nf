/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


/*
*******************************************************************************
* version_check
*******************************************************************************
 */
def version_check(required_ver, current_ver){

	try {
	    if( ! current_ver.matches(">= $required_ver") ){
		throw GroovyException('Nextflow version too old')
	    }
	} catch (all) {
		log.error "Error: Nextflow version $required_ver required! " +
		      "You are running version $current_ver.\n" +
		      "Please update Nextflow.\n"
		exit 1

	}
}


/*
*******************************************************************************
* Function to show yes/no prompt
*******************************************************************************
 */
def prompt(input){

	if(input == "n"){
		exit 1
	}
	if(input == "y"){
		return(true)
	}
	if(input != "n" || input!= "y"){
		println "Please use 'y' for yes and 'n' for no."
		prompt(System.console().readLine 'Do you want to continue again (y/n)?')
	}
}


/*
*******************************************************************************
* Nextpie
*******************************************************************************
 */

@Grab('io.github.http-builder-ng:http-builder-ng-okhttp:0.14.2')
@Grab(group='org.slf4j', module='slf4j-api', version='1.7.32')

import static groovy.json.JsonOutput.toJson
import static groovyx.net.http.HttpBuilder.configure

import groovyx.net.http.*
import static groovyx.net.http.MultipartContent.multipart

def Nextpie(host, port, traceFile, Workflow, Version, Group, Project, APIkey){
    File myFile = new File(traceFile)

    try {
	    def posts = configure {
		request.uri = 'http://'+host+':'+port
		request.uri.path = '/api/v1.0/upload-data'
		request.headers['X-API-KEY'] = APIkey
		request.contentType = 'multipart/form-data'
		request.body = multipart {
		    field 'Workflow', Workflow
		    field 'Version', Version
		    field 'Group', Group
		    field 'Project', Project
		    part 'File', 'Trace.txt', 'text/plain', myFile
		}
		request.encoder 'multipart/form-data', OkHttpEncoders.&multipart

		response.success { FromServer fs, Object body ->
		    return  body
		}

		response.when(401){ FromServer fs, Object body ->
			return "UNAUTHORIZED (401)"
		}

		response.when(403){ FromServer fs, Object body ->
			return "FORBIDDEN (403)"
		}
		response.when(404){ FromServer fs, Object body ->
			return "NOT FOUND (404)"
		}
	    }.post()
    }
    catch (Exception ce){
    	return ce
    }
}
