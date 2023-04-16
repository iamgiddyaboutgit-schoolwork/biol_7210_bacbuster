#Insert Shebang?

#Might need to use virtualenv if conda doesn't cooperate

#Generate flask template for later builds 
from flask import Flask, render_template, url_for

app = Flask(__name__)
@app.route('/')

def index():
    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)

