import streamlit as st

#Provide a large title in bold.
st.text("#BacBuster Isolate Analyzer")

#Provide some info on how the app works and the authors.
st.text("##This app seeks to simplify the process of ")

#Create form so that uploads and email are all processed together.
with st.form("func"):

    #Give the user the ability to upload isolates TODO Increase file size limit and theming.
    upload = st.file_uploader("Submit raw paired reads of the isolate(s).", accept_multiple_files=True)

    #Find out how many pairs of files there are.
    pairCount = len(upload) / 2

    #Add a radio button that allows the user to decide on how much of the pipeline should be run.
    st.radio("Desired Output", ["Genome Assembly", "Gene Prediction", "Functional Annotation"])

    #Add a text input for email. TODO make this optional.
    email = st.text_input("Provide email for to receive results later.")

    #Submit button for the form.
    submitted = st.form_submit_button()
    if submitted:
        st.write("Inputs provided. Please wait for results!\n Note that pressing the submit button again will reset the pipeline.")

#RUN THE NEXTFLOW PIPELINE AND RETURN RESULTS
import subprocess

#Assembly Step


#Notify user that assembly is finished.


#Prediction Step


#Notify user that prediction is finished.
#Prediction step completed 2/3, hold tight!

#Annotation Step

#Write to results screen. Before completion, have some sort of progress printed to user.

#Provide download button for the results.
#db = st.download_button()

#Email results to input email using smtp. Wait until the app has been tested without.

