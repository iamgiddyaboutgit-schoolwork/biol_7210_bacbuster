#!/usr/bin/env python3

from flask import Flask
from flask import request

app = Flask(__name__)

# This is a view function that tells the Flask app instance
# what code to run.  Navigation on the website to the route
# causes this view function to execute.  This view function
# returns a response to the client. 
# Through the magic of contexts, view functions like the one
# below can get access to global objects.  Contexts are nice
# because they don't cause clashes across threads and multiple
# users.
@app.route("/")
def index():
    user_agent = request.headers.get("User-Agent")
    return f"<h1>BacBuster will be used on your browser of {user_agent}.</h1>"


#####################################################################
if __name__ == '__main__': 
    # Do not use run() in a production setting. It is not intended 
    # to meet security and performance requirements for a production 
    # server.
    app.run(debug=True)