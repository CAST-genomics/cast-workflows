set -e
DX_COMPILER_JAR="/Users/yang/Desktop/CAST_repository/WDL_tools/dxCompiler-2.10.4.jar"
PROJID="project-GbjB4x0JX3JyK5PZK67q0fZG"
# java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/trim_vcf.wdl
# java -jar "$DX_COMPILER_JAR" compile ../wdl/trim_vcf.wdl -project $PROJID -folder /yal084/applets/trim_vcf -streamFiles all -archive
java -jar ~/Desktop/CAST_repository/WDL_tools/wdlTools-0.17.12.jar check ../wdl/concat_vcf.wdl
java -jar "$DX_COMPILER_JAR" compile ../wdl/concat_vcf.wdl -project $PROJID -folder /yal084/applets/concat_vcf -streamFiles all -archive


