
# coding: utf-8

# In[ ]:

def sendtext(scriptfile, outputfile, from_addr="ptexter4242@gmail.com", phonenumber="4259620542", carriername='verizon', login='ptexter4242@gmail.com', password='2424RETXETp', smtpserver='smtp.gmail.com:587'):

    import smtplib

    def carrier(carriername):
        carriers = {'alltel': '@message.alltel.com',
                    'att': '@txt.att.net',
                    'boost': '@myboostmobile.com',
                    'cingular': '@cingularme.com',
                    'cingular2': '@mobile.mycingular.com',
                    'nextel': '@messaging.nextel.com',
                    'sprint': '@messaging.sprintpcs.com',
                    'tmobile': '@tmomail.net',
                    'verizon': '@vtext.com',
                    'virgin': '@vmobl.com'}

        return carriers[carriername]



    to_email = phonenumber + carrier(carriername)

    server = smtplib.SMTP(smtpserver)
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login(login,password)

    problems = server.sendmail(from_addr, to_email, "Howssss it? " + scriptfile + " has finished and the output is in " + outputfile)

    server.quit()

