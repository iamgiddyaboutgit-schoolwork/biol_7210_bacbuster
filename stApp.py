import streamlit as st

#Provide a large title in bold.
st.text("BacBuster Isolate Analyzer")

#Provide some info on how the app works and the authors.


#Create form so that uploads and email are all processed together.
with st.form("func"):

    #Give the user the ability to upload isolates TODO Increase file size limit.
    upload = st.file_uploader("Submit raw paired reads of the isolate(s).", accept_multiple_files=True)

    #Find out how many pairs of files there are.
    pairCount = len(upload) / 2

    #Add a text input for email. TODO make this optional.
    email = st.text_input("Provide email for to receive results later.")

    #Submit button for the form.
    submitted = st.form_submit_button
    if submitted:
        st.write("Inputs provided. Please wait for results!")



#Provide download button for the results.
#db = st.download_button()

#Email results to input email using smtp.
