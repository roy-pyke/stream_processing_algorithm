#!/usr/bin/env bash
# Fetch and build the journal version's code (bounded deletions), unmodified, plus our I/O harness.
#   ext/omni-deletes   github.com/wiegerpunter/omni-deletes at commit 17b969e
#   ext/jars           the dependencies of their pom.xml that the PVLDB code does not ship: Guava 33.3.1-jre,
#                      failureaccess 1.0.2, Jackson 2.18.0 (Maven Central, SHA-1 checked)
#   ext/build/deletes  class files (javac -sourcepath: only the classes our harness reaches)
# Needs scripts/setup_upstream.sh first (JDK in ext/java_home.txt, dependency jars of the PVLDB code).
set -euo pipefail
cd "$(dirname "$0")/.."
COMMIT=17b969e
[ -d ext/omni-deletes/.git ] || git clone --quiet https://github.com/wiegerpunter/omni-deletes ext/omni-deletes
git -C ext/omni-deletes checkout --quiet "$COMMIT"
mkdir -p ext/jars ext/build/deletes
fetch() {
  local f=ext/jars/$(basename "$1")
  [ -s "$f" ] || curl -sSfL -o "$f" "https://repo1.maven.org/maven2/$1"
  echo "$2  $f" | shasum -a 1 -c - > /dev/null || { echo "checksum mismatch: $f"; exit 1; }
}
fetch com/google/guava/guava/33.3.1-jre/guava-33.3.1-jre.jar 852f8b363da0111e819460021ca693cacca3e8db
fetch com/google/guava/failureaccess/1.0.2/failureaccess-1.0.2.jar c4a06a64e650562f30b7bf9aaec1bfed43aca12b
fetch com/fasterxml/jackson/core/jackson-databind/2.18.0/jackson-databind-2.18.0.jar 8dba1f789a75fc30b59303574fe2b269afa4d3bc
fetch com/fasterxml/jackson/core/jackson-core/2.18.0/jackson-core-2.18.0.jar 65e8ead7de5d8f7a53e296c363bea3182f21f925
fetch com/fasterxml/jackson/core/jackson-annotations/2.18.0/jackson-annotations-2.18.0.jar 9bddcc56af9d90f902ef5dba7348102cd12b04e2
JH=$(cat ext/java_home.txt)
CP="$(ls ext/jars/*.jar | tr '\n' ':')$(ls ext/upstream/out/artifacts/DSCM_jar/*.jar | grep -v '/DSCM\.jar$' | tr '\n' ':')"
"$JH/bin/javac" -nowarn -d ext/build/deletes -cp "$CP" -sourcepath ext/omni-deletes/src/main/java repro/omni_deletes/DeletesHarness.java
echo "ext/build/deletes:$CP" > ext/build/deletes_cp.txt
echo "built: java -cp \$(cat ext/build/deletes_cp.txt) omni.DeletesHarness ..."
