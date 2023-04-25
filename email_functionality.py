import os
import requests
import streamlit as st

# Function to generate output file

# Function to send email
def send_email(api_key, domain, recipient_email, subject, body, files=[]):
    # Construct the Mailgun API endpoint
    endpoint = f"https://api.mailgun.net/v3/{domain}/messages"

    # Create a dictionary of the email message data
    message_data = {
        "from": f"Streamlit <mailgun@{domain}>",
        "to": recipient_email,
        "subject": subject,
        "text": body,
    }

    # Add file attachments
    files = [("attachment", (os.path.basename(file), open(file, "rb"))) for file in files]

    # Send the email using the Mailgun API
    response = requests.post(endpoint, auth=("api", api_key), files=files, data=message_data)

    # Check if the email was successfully sent
    if response.status_code == 200:
        st.success("Email sent successfully!")
    else:
        st.error(f"Failed to send email. Status code: {response.status_code}. Response: {response.content}")

# Streamlit app
def app():
    st.title('Email Output File')

    # Get recipient email
    recipient_email = st.text_input('Enter recipient email')

    # Generate output file
    output_file_path = "/home/sandhya/output_file.txt"

    # Send email
    if st.button('Send email'):
        # Get Mailgun API key and domain from environment variables
        api_key = "d549f8b9c018df234c32d40f87d87fa8"
        print(api_key)
        domain = "sandbox6108f1c320074d4ca0df5f7a5054c05c.mailgun.org"
        print(domain)
        # Send email
        send_email(api_key, domain, recipient_email, 'Output File', 'Here is your output file', files=[output_file_path])

# Run the app
if __name__ == '__main__':
    app()

