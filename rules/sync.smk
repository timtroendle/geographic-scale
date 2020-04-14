# Rules to sync to and from Euler

EULER_URL = "euler.ethz.ch"
EULER_BASE_DIR = "~/Develop/geographical-scale/"
EULER_BUILD_DIR = EULER_BASE_DIR + "build/"
LOCAL_EULER_RESULTS = "./build/euler/"
LOCAL_PUBLICATION_RESULTS = "./build/euler/publish/"
SYNCIGNORE = ".syncignore"
SYNCIGNORE_BUILD = ".syncignore-build"
SYNCIGNORE_PUBLISH = ".syncignore-publish"


rule send:
    message: "Send changes to Euler"
    shell:
        "rsync -avzh --progress --delete -r . --exclude-from={SYNCIGNORE} {EULER_URL}:{EULER_BASE_DIR}"


rule receive:
    message: "Receive build changes from Euler"
    shell:
        "rsync -avzh --progress --delete -r --exclude-from={SYNCIGNORE_BUILD} {EULER_URL}:{EULER_BUILD_DIR} {LOCAL_EULER_RESULTS}"


rule prepare_publication:
    message: "Receive build changes from Euler as preparation for publication"
    shell:
        "rsync -avzhH --progress --delete -r --exclude-from={SYNCIGNORE_PUBLISH} {EULER_URL}:{EULER_BUILD_DIR} {LOCAL_PUBLICATION_RESULTS}"


rule clean_euler:
    message: "Clean results downloaded from Euler"
    shell:
        """
        rm -r {LOCAL_EULER_RESULTS}/*
        """
