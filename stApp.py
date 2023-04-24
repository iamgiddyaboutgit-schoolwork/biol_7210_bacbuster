import streamlit as st
import math
#import sys
import subprocess

#Provide a large title in bold.
st.title(":boom: BacBuster Isolate Analyzer :boom:")

#Provide some info on how the app works and the authors.
st.text("This app simplifies analysis of raw paired reads from outbreak isolates.\nUsers can process the samples up to the AMR annotation, but can stop before.")


#Create form so that uploads and email are all processed together.
with st.form("func"):

    #Give the user the ability to upload isolates
    upload = st.file_uploader("Upload raw paired reads of the isolate(s).", accept_multiple_files=True)


    #Add a radio button that allows the user to decide on how much of the pipeline should be run.
    stop = st.radio("Desired Output", ["1. Genome Assembly", "2. Gene Prediction", "3. Functional Annotation"])
    st.write("⚠️ NOTE: Selecting a later step will increase the runtime and size of the final results folder.")

    #Add a text input for email. TODO make this optional.
    email = st.text_input("Provide email if you would like to receive results later.", placeholder="kingjordan@gmail.com")

    #Provide the user with additional parameters to influence the assembly/prediction/annotation steps?


    #Submit button for the form. TODO disable on click? A mistaken submission would need to reopen the app then.
    submitted = st.form_submit_button()
    if submitted:
        st.write("Inputs provided. Please wait for results!\n Note that pressing the submit button again will reset the pipeline.")



#RUN THE NEXTFLOW PIPELINE AND RETURN RESULTS USING SUBPROCESS.
#import subprocess


files = subprocess.run(["nextflow run", "src/nextflow/BacBuster.nf"])
#files = subprocess.run([f"{sys.executable}", "src/nextflow/BacBuster.nf"])
st.write(files.stdout)

#Generate the folder that nextflow will read from. Let the whole thing run for now, but later introduce break points based on choice.
#subprocess.run("mkdir")

#THESE WIDGETS SHOULD ONLY BE DRAWN AFTER SUBMIT IS PRESSED.
#Find out how many pairs of files there are.
pairCount = math.floor(len(upload) / 2) # if no files, returns 0.

#Check that this is at least 1 pair (the user has input at least 2 files.)
if submitted and pairCount == 0:
    st.error("2 files are needed at minimum, please confirm that at least one pair of reads has been uploaded!", icon="🚨")


#Assembly Step
#subprocess(nextflow run )

assemblyCheck = False #Flip this flag once Nextflow pipeliine is done AND results from this section are not empty.

#Notify user that assembly is finished.
if submitted and assemblyCheck:
    st.success(' Assembly completed!', icon="✅")

#Show any results that can be visualized from this step.

#Prediction Step


#Notify user that prediction is finished.

#Show any results that can be visualized from this step.

predictionCheck = False #Flip this flag once Nextflow pipeliine is done.
if submitted and predictionCheck:
    st.success(' Prediction completed! If Functional Annotation was selected, hold tight! Downloadableesults incoming.', icon="✅")
#Prediction step completed 2/3, hold tight!

#Annotation Step

#Write to results screen. Before completion, have some sort of progress printed to user.


#Show any results that can be visualized from this step.

#st.write() handles many data inputs, such as images.


#Provide download button for the results. TODO Hide this until the pipeline has finished running later.
#db = st.download_button("Download zipped final results.", f, file_name="results.zip")

#Email results to input email using smtp. Wait until the app has been tested without.
if submitted and len(email) > 0:
    #Send email
    st.write("email provided")


#TODO: Add expander with more information on the workflow + flowchart.
#Another TODO: Use session states over a check for submitted for showing the results.

st.text("BacBuster was developed during completion of BIOL 7210, Computational Genomics, at Georgia Tech.")