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

import "tasks/utilities.wdl" as utilities
import "tasks/minimap2.wdl" as minimap2
import "tasks/sambamba.wdl" as sambamba
import "tasks/clair.wdl" as clair
import "tasks/GATK4.wdl" as GATK4
import "tasks/bcftools.wdl" as bcftools
import "tasks/bgzip.wdl" as bgzip
import "tasks/tabix.wdl" as tabix
import "tasks/seqkit.wdl" as seqkit

workflow ONT_VariantCalling {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.4"
		date: "2021-05-10"
	}

	input {
		String fastqPath
		String outputPath

		File refFa
		File refFai

		String extFq = "*.fastq"
		Int maxDepth = 1
		Boolean gzip = false

		String modelPath

		Int qual = 748
		File? bedRegions

		Int? maxLen
		Int? maxQual
		Int? minLen
		Int? minQual

		String? filterRegionVCF
		String? filterIncludeVCF

		String? name
	}

	String sampleName = if defined(name) then "~{name}" else "sample"


################################################################################
## Alignment

	call utilities.findFiles as FINDFILES {
		input :
			path = fastqPath,
			regexpName = extFq,
			maxDepth = maxDepth
	}

	# seqkt filter

	String extConc = if gzip then ".fastq.gz" else ".fastq"

	call utilities.concatenateFiles as CONCATENATEFILES {
		input :
			in = FINDFILES.files,
			name = sampleName + extConc,
			gzip = gzip,
			outputPath = outputPath + "/fastq_concatenate/"
	}

	if (defined(maxLen) || defined(maxQual) || defined(minLen) || defined(minQual)) {
		call seqkit.seq_filter as FQ_FILTER {
			input :
				in = CONCATENATEFILES.outputFile,
				maxLen = maxLen,
				maxQual = maxQual,
				minLen = minLen,
				minQual = minQual,
				subString = extConc,
				subStringReplace = ".filtered" + extConc,
				outputPath = outputPath + "/fastq_concatenate/"
		}
	}

	call minimap2.mapOnt as MAPONT {
		input :
			fastq = select_first([FQ_FILTER.outputFile, CONCATENATEFILES.outputFile]),
			refFasta = refFa,
			sample = sampleName,
			outputPath = outputPath + "/Alignment/"
	}

	call sambamba.index as INDEX {
		input :
			in = MAPONT.outputFile,
			outputPath = outputPath + "/Alignment/"
	}

################################################################################

################################################################################
## Variant Calling

	call utilities.fai2bed as FAI2BED {
		input :
			in = refFai
	}

	call utilities.bed2Array as BED2ARRAY {
		input :
			bed = select_first([bedRegions,FAI2BED.outputFile])
	}

	scatter (region in BED2ARRAY.bedObj) {
		call clair.callVarBam as CALLVARBAM {
			input :
				modelPath = modelPath,
				refGenome = refFa,
				refGenomeIndex = refFai,
				bamFile = MAPONT.outputFile,
				bamFileIndex = INDEX.outputFile,
				name = "~{sampleName}-~{region['chrom']}_~{region['start']}_~{region['end']}",
				sampleName = sampleName,
				qual = qual,
				outputPath = "~{outputPath}/Variant-Calling/clair/scatter/",
				contigName = region["chrom"],
				ctgStart = region["start"],
				ctgEnd = region["end"],
				delay = 10
		}
	}

	call GATK4.gatherVcfFiles as GATHERVCFFILES {
		input :
			in = CALLVARBAM.outputFile,
			outputPath = outputPath + "/Variant-Calling/clair/",
			subString = "-(chr)?[0-9WXYMT]+_[0-9]+_[0-9]+\.clair\.vcf$"
	}

	# filter
	# gz and index
	call bgzip.compress as BGZ_VCF {
		input :
			in = GATHERVCFFILES.outputFile,
			outputPath = outputPath + "/Variant-Calling/clair/"
	}
	call tabix.index as TBI {
		input :
			in = BGZ_VCF.outputFile,
			outputPath = outputPath + "/Variant-Calling/clair/"
	}

	if (defined(filterRegionVCF) || defined(filterIncludeVCF)) {
		call bcftools.view as VIEW {
			input :
				in = BGZ_VCF.outputFile,
				index = TBI.outputFile,
				region = filterRegionVCF,
				includeExp = filterIncludeVCF,
				subStringReplace = ".filtered",
				outputPath = outputPath + "/Variant-Calling/clair/"
		}
	}

################################################################################

################################################################################
## Outputs

	output {
		File bam = MAPONT.outputFile
		File bai = INDEX.outputFile
		File vcf = GATHERVCFFILES.outputFile
		File vcfIdx = GATHERVCFFILES.outputFileIdx
	}

################################################################################
}
