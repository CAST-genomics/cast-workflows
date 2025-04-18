set -e
DX_COMPILER_JAR="/Users/yang/Desktop/CAST_repository/WDL_tools/dxCompiler-2.10.4.jar"
PROJID="project-GbjB4x0JX3JyK5PZK67q0fZG"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/trim_vcf.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/trim_vcf.wdl -project $PROJID -folder /yal084/applets/trim_vcf -streamFiles all -archive
#
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/concat_vcf.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/concat_vcf.wdl -project $PROJID -folder /yal084/applets/concat_vcf -streamFiles all -archive

# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/bfileTovcf.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/bfileTovcf.wdl -project $PROJID -folder /yal084/applets/bfileTovcf -streamFiles all -archive

# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/bgenTovcf.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/bgenTovcf.wdl -project $PROJID -folder /yal084/applets/bgenTovcf -streamFiles all -archive

# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/hg19Tohg38.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/hg19Tohg38.wdl -project $PROJID -folder /yal084/applets/hg19Tohg38 -streamFiles all -archive

# for file in $(ls -d ../wdl/*);do 
#   echo $(basename $file)
# done

# echo "bfileToVCF workflow-id"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/bfileTovcf.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/bfileTovcf.wdl -project $PROJID -folder /yal084/applets/bfileTovcf -streamFiles all -archive

# echo "hg19Tohg38 workflow-id"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/hg19Tohg38.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/hg19Tohg38.wdl -project $PROJID -folder /yal084/applets/hg19Tohg38 -streamFiles all -archive

# echo "splitVCF workflow-id"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/splitVCF.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/splitVCF.wdl -project $PROJID -folder /yal084/applets/splitVCF -streamFiles all -archive
#
# echo "imputation workflow-id"
# java -jar ~/desktop/cast_repository/wdl_tools/wdltools-0.17.12.jar check ../wdl/batch_imputation_bref3.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/batch_imputation_bref3.wdl -project $PROJID -folder /yal084/applets/batch_imputation_bref3 -streamfiles all -archive


# echo "extract_trs workflow"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/extract_trs.wdl
# echo "check passed, start compiling"
# java -jar "$DX_COMPILER_JAR" compile ../wdl/extract_trs.wdl -project $PROJID -folder /yal084/applets/extract_trs -streamFiles all -archive

echo "merge_and_annotation workflow"
java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/merge_and_annotation.wdl
echo "check passed, start compiling"
java -jar "$DX_COMPILER_JAR" compile ../wdl/merge_and_annotation.wdl -project $PROJID -folder /yal084/applets/merge_and_annotation -streamFiles all -archive

# echo "bcf_liftover workflow-id"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/bcf_liftover.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/bcf_liftover.wdl -project $PROJID -folder /yal084/applets/bcf_liftover -streamFiles all -archive

# echo "splitVCF_para workflow-id"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/splitVCF_para.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/splitVCF_para.wdl -project $PROJID -folder /yal084/applets/splitVCF_para -streamFiles all -archive

