import streamlit as st

#Provide a large title in bold.
st.title(":boom: BacBuster Isolate Analyzer :boom:")

#Provide some info on how the app works and the authors.
st.text("This app seeks to simplify the process of analyzing raw paired reads from outbreak isolates. Users can ask for a readout of antimicrobial genes identified, or stop the pipeline at a preliminary step.")

#Create form so that uploads and email are all processed together.
with st.form("func"):

    #Give the user the ability to upload isolates TODO Increase file size limit and theming.
    upload = st.file_uploader("Submit raw paired reads of the isolate(s).", accept_multiple_files=True)


    #Add a radio button that allows the user to decide on how much of the pipeline should be run.
    st.radio("Desired Output", ["1. Genome Assembly", "2. Gene Prediction", "3. Functional Annotation"])
    st.write("⚠️ NOTE: Selecting a later step will increase the runtime and size of the final results folder.")

    #Add a text input for email. TODO make this optional.
    email = st.text_input("Provide email for to receive results later.")

    #Submit button for the form.
    submitted = st.form_submit_button()
    if submitted:
        st.write("Inputs provided. Please wait for results!\n Note that pressing the submit button again will reset the pipeline.")

#RUN THE NEXTFLOW PIPELINE AND RETURN RESULTS
import subprocess

#Find out how many pairs of files there are.
pairCount = len(upload) / 2 # if no files, returns 0.
#Check that this is at least 1 (the user has input at least 2 files.)
if pairCount == 0:
    st.error("2 files are needed at minimum, please confirm that at least one pair of reads has been uploaded!", icon="🚨")


#Assembly Step


#Notify user that assembly is finished.
st.success('Assembly completed!', icon="✅")

#Prediction Step


#Notify user that prediction is finished.

st.success('Prediction completed! If Functional Annotation was selected, hold tight! Results incoming.', icon="✅")
#Prediction step completed 2/3, hold tight!

#Annotation Step

#Write to results screen. Before completion, have some sort of progress printed to user.

#Provide download button for the results.
#db = st.download_button()

#Email results to input email using smtp. Wait until the app has been tested without.
if email is not None:
    #Send email
    print("email provided")