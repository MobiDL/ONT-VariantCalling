version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2021 MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import "tasks/utilities.wdl" as runUtilities
import "tasks/fastqc.wdl" as runFastqc
import "tasks/GATK4.wdl" as runGATK
import "tasks/multiqc.wdl"as runMultiqc
import "tasks/minimap2.wdl" as runMinimap2
import "tasks/sambamba.wdl" as runSambamba
import "tasks/clair3.wdl" as runClair3
import "tasks/nanorepeat.wdl" as runNanorepeat
import "tasks/bcftools.wdl" as runBcftools

workflow ONT_2023 {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.1.0"
		date: "2023-10-18"
	}

	input {
		String fastqPath
		String outputPath

		File refFa
		File refFai
		File refDict

		String extFq = "*.fastq"
		Int maxDepth = 1
		Boolean gzip = false

		String modelPath
		File bedRegions
		File nanorepeatBedRegions = bedRegions
		
		String? nanorepeatDataType

		String? name
		Int maxThreads = 1
	}

	String sampleName = if defined(name) then "~{name}" else "sample"


################################################################################
## Fastq

	call runUtilities.findFiles as FINDFILES {
		input :
			path = fastqPath,
			regexpName = extFq,
			maxDepth = maxDepth
	}

	String extConc = if gzip then ".fastq.gz" else ".fastq"

	call runUtilities.concatenateFiles as CONCATENATEFILES {
		input :
			in = FINDFILES.files,
			name = sampleName + extConc,
			gzip = gzip,
			outputPath = outputPath + "/fastq/"
	}

	call runFastqc.fastqc as FASTQC {
		input:
			in = [CONCATENATEFILES.outputFile],
			threads = maxThreads,
			outputPath = outputPath+ "/qc/"
	}


################################################################################
## Bed

	# call runGATK.bedToIntervalList as BED2INTERVAL {
	# 	input:
	# 		in = bedRegions,
	# 		refDict = refDict,
	# 		outputPath = outputPath+ "/intervals/"
	# }

################################################################################
## Alignement

	call runMinimap2.mapOnt as MAPONT {
		input :
			fastq = CONCATENATEFILES.outputFile,
			refFasta = refFa,
			sample = sampleName,
			threads = maxThreads,
			outputPath = outputPath + "/alignment/"
	}

	call runSambamba.index as INDEX {
		input :
			in = MAPONT.outputFile,
			threads = maxThreads,
			outputPath = outputPath + "/alignment/"
	}

################################################################################
## QC

	call runGATK.depthOfCoverage as DEPTHOFCOVERAGE {
		input:
			in = MAPONT.outputFile,
			bamIndex = INDEX.outputFile,
			referenceFasta = refFa,
			referenceFai = refFai,
			referenceDict = refDict,
			outputPath = outputPath + "/qc/",			
			intervals = bedRegions
	}


################################################################################
## Variant Calling
	call runClair3.clair3 as CLAIR3 {
		input :
			modelPath = modelPath,
			refGenome = refFa,
			refGenomeIndex = refFai,
			bamFile = MAPONT.outputFile,
			bamFileIndex = INDEX.outputFile,
			bedFile = bedRegions,
			threads = maxThreads,
			outputPath = outputPath + "/clair3/"
	}

	if (defined(nanorepeatDataType)) {
		call runNanorepeat.nanorepeat as NANOREPEAT {
			input:
				bamFile = MAPONT.outputFile,
				bamFileIndex = INDEX.outputFile,
				refGenome = refFa,
				refGenomeIndex = refFai,
				bedFile = nanorepeatBedRegions,
				dataType = nanorepeatDataType,
				threads = maxThreads * 2,
				outputPath = outputPath + "/nanorepeat/" + sampleName
		}
	}


################################################################################
### normalise VCF

	call runBcftools.norm as BCFTOOLSSPLIT {
		input:
			in = CLAIR3.outputFile,
			outputPath = outputPath,
			splitMA = true,
			multiallelicType = "both",
			refFasta = refFa
        }

################################################################################
## MultiQC

	call runMultiqc.multiqc as MULTIQC {
		input:
			MetrixFiles = FASTQC.outHTML,
			outputPath = outputPath,
			path_to_check = outputPath
	}

################################################################################
## Outputs

	output {
		# Array[File] qc = FASTQC.outHTML
		File bam = MAPONT.outputFile
		File bai = INDEX.outputFile
		File? vcf = CLAIR3.outputFile
		File multiqc = MULTIQC.outputFile
	}

################################################################################
}
