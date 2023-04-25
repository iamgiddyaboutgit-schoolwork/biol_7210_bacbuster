import streamlit as st
import math
import os
import subprocess

#DRAW FRONT-END
#Provide a large title in bold.
st.title(":boom: BacBuster Isolate Analyzer :boom:")

#Provide some info on how the app works and the authors.
st.text("This app simplifies analysis of raw paired reads from outbreak isolates.\nUsers can process the samples up to the AMR annotation.")

#Create form so that uploads and email are all processed together.
with st.form("func"):

    #Give the user the ability to upload isolates
    upload = st.file_uploader("Upload raw paired reads of the isolate(s) as _fq.gz", accept_multiple_files=True)

    #Store the isolates in a new directory: inDir that can be passed later.
    #if len(inDir)
    #    for read in upload:
    #        with open(read.name, "wb") as f:
    #            f.write()

    #Add a radio button that allows the user to decide on how much of the pipeline should be run.
    #stop = st.radio("Desired Output", ["1. Genome Assembly", "2. Gene Prediction", "3. Functional Annotation"])
    #st.write("⚠️ NOTE: Selecting a later step will increase the runtime and size of the final results folder.")

    #Add a text input for email.
    email = st.text_input("Provide email if you would like to receive results later.", placeholder="kingjordan@gmail.com")

    #Provide the user with additional parameters to influence the assembly/prediction/annotation steps?


    #Submit button for the form.
    submitted = st.form_submit_button()
    if submitted:
        st.write("Please wait for results!\n Note that pressing the submit button again will reset the pipeline.")

#EVALUATE USER INPUT
#Find out how many pairs of files there are.
#pairCount = math.floor(len(upload) / 2) # if no files, returns 0.
pairCount = len(upload) // 2

#Check that this is at least 1 pair (the user has input at least 2 files.)
if submitted and pairCount == 0:
    st.error("Two files are needed at minimum, please confirm that at least one pair of reads has been uploaded! Then press Submit again.", icon="🚨")
    submitted = False

#Check that an even number of reads has been provided.
if submitted and len(upload) % 2 == 1:
    st.error("Odd number of files detected, please confirm that each isolate is submitted as a pair. Then press Submit again.", icon="🚨")
    submitted = False

#RUN THE NEXTFLOW PIPELINE AND RETURN RESULTS USING SUBPROCESS.

#files = subprocess.run(["nextflow", "run", "BacBuster.nf", "--seq_reads", "/home/andy/compGen/Team3-WebServer/testing_data/sequencing_reads"], capture_output=True, text=True)
#TODO Implement user input.
if submitted: #All the steps that occur once the user input has been verified.
    files = subprocess.run(["nextflow", "run", "BacBuster.nf", "--seq_reads" "testing_data/sequencing_reads"], capture_output=True, text=True)
    #files = subprocess.run(["nextflow", "run", "BacBuster.nf", "--seq_reads", "inDir"], capture_output=True, text=True)
#files = subprocess.run([f"{sys.executable}", "src/nextflow/BacBuster.nf"])
    st.write(files.stdout)

#Generate the folder that nextflow will read from. Let the whole thing run for now, but later introduce break points based on choice.
#subprocess.run("mkdir")


#assemblyCheck = False #Flip this flag once Nextflow pipeliine is done AND results from this section are not empty.

#Notify user that assembly is finished.
#if submitted and assemblyCheck:
#    st.success(' Assembly completed!', icon="✅")

#Show any results, if possible.
#Show any results that can be visualized from prediction.
#predictionCheck = False #Flip this flag once Nextflow pipeliine is done.
#if submitted and predictionCheck:
#    st.success(' Prediction completed! If Functional Annotation was selected, hold tight! Downloadable results incoming.', icon="✅")

#Write to results screen. Before completion, have some sort of progress printed to user.

#st.write() handles many data inputs, such as images.


#Provide download button for the results. TODO Hide this until the pipeline has finished running later.

#Generate zip file of results.
    #subprocess.run(["rm", "-r", "work/conda"])
    #subprocess.run(["zip", "-r", "temp.zip", "work"])
    subprocess.run(["zip", "-r", "temp.zip", "output"])
    st.success("Pipeline and packaging finished! Results below", icon="✅")

#Draw a button and pass the zip.
    with open("temp.zip", "rb") as fp:
        db = st.download_button(label="Download zipped final results.", data = fp, file_name="results.zip", mime="application/zip")

#Email results to input email using smtp. Wait until the app has been tested without.
    if len(email) > 0:
        #Send email
        st.write("Email provided!")

    #Add a bunch of imports for the email package.


#TODO: Add expander with more information on the workflow + flowchart.
#Another TODO: Use session states over a check for submitted for showing the results.

st.text("BacBuster was developed for BIOL 7210, Computational Genomics, at Georgia Tech\nduring the Spring 2023 semester. Thanks for trying it!")
