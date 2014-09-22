main() {
    dx-download-all-inputs

    # Download with the streaming utility, and compare
    mkdir -p tmp
    dx-stream-input seq1 > tmp/seq1
    DIFF=$(diff tmp/seq1 in/seq1/*)
    if [ "$DIFF" != "" ]
    then
        echo "Mismatch between streaming and download-all"
        exit 1
    fi

    mkdir -p out/foo
    echo "ABCD" > out/foo/A.txt

    dx-upload-all-outputs

    # use streaming to compare the two upload options
}
