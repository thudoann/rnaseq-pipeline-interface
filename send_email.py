import smtplib
from email.mime.text import MIMEText
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.DEBUG)


def send_email(subject, body, to_email):
    from_email = 'rnaseq.pipeline@gmail.com'
    smtp_server = 'smtp.gmail.com'
    smtp_user = 'rnaseq.pipeline@gmail.com'
    smtp_port = 587
    smtp_password = 'lwyhliurwyclqxpn'

    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = from_email
    msg['To'] = to_email

    logging.debug(f"Email details: Subject: {subject}, From: {from_email}, To: {to_email}")
    logging.debug(f"SMTP Server: {smtp_server}, Port: {smtp_port}")

    try:
        with smtplib.SMTP(smtp_server, smtp_port) as server:
            server.starttls()
            server.login(smtp_user, smtp_password)
            server.sendmail(from_email, to_email, msg.as_string())
        logging.info(f"Email sent successfully to {to_email}")
    except Exception as e:
        logging.error(f"Failed to send email: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python send_email.py <summary_file> <recipient_email>")
        sys.exit(1)

    summary_file = sys.argv[1]
    recipient_email = sys.argv[2]

    logging.debug(f"Summary file: {summary_file}, Recipient email: {recipient_email}")

    try:
        with open(summary_file, 'r') as f:
            summary_content = f.read()
        logging.debug(f"Summary content:\n{summary_content}")
        send_email("My pipeline execution", summary_content, recipient_email)
    except Exception as e:
        logging.error(f"Failed to read summary file: {e}")