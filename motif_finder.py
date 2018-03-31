#!/usr/bin/env python

from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO
import re
from tkinter import *

#this is the regex string used to find the motif 
motif_regex = r"[^P][^PKRHW][VLSWFNQ][ILTYWFN][FIY][^PKRH]"
#this is the regex string used to validate the input so it only allows valid swissprot accession numbers
accession_regex = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"

#the found motif sequences, end and start positions will be placed in this dictionary as key:value pairs
sequences_and_positions = dict()

#the user interface of this program is the tkinter UI library which comes standard with Python
#here we create the window and give it a title
window = Tk()
window.wm_title("motif locator")

#this update function will be run every time the "Find motifs" button is clicked
def update():
    #get the protein code from the text entry box
    protein = protein_number.get()
    if re.search(accession_regex,protein) is None:
        #if the regex string returns nothing, we need to tell the user so they can enter a new one
        #we must also end the function so it doesn't try to parse a null object
        #to alert the user, we change the text in the entry box to an alert message
        protein_number.delete(0, END)
        protein_number.insert(0, "Not a valid protein number")
        return
    #if the protein number is valid, the string is requested from Swissprot and parsed
    protein = protein_number.get()
    handle = ExPASy.get_sprot_raw(protein)
    record = SwissProt.read(handle)
    handle.close()
    description = record.description
    sequence = record.sequence
    runs = re.finditer(motif_regex,sequence)
    #clear the dictionary so old values arent placed in the list of found sequences
    sequences_and_positions = dict()
    for match in runs:
        """in this loop, we go through the list of found sequences, and put the sequence into
           the sequences_and_positions dictionary as the key value, its corresponding value is a tuple
           which contains the start and end position of the sequence in the protein
           when the UI updates, another loop will use this dictionary to update the console
        """
        run_start = match.start()
        run_end = match.end()
        #the sequence is between the start and end positions
        seq = sequence[run_start:run_end]
        sequences_and_positions[seq] = (run_start+1,run_end+1)
    """in the following lines the UI is updated
       first all content in the console (a Tkinter listbox object) is deleted
       then a "header" is added which indicates that the 3 columns are the sequence, start and end
       respectively, then a loop goes through every key:value pair and puts the key, and two values
       in the corresponding tuple in the listbox. Each value has a padding of 8 so the items line up.
    """
    info_frame.delete(1.0,END)
    info_frame.insert(END,"Sequence    Start   End\n")
    for key in sequences_and_positions:
        begin = sequences_and_positions[key][0]
        end = sequences_and_positions[key][1]
        #the string format function has a padding option, {:8} means each item takes 8 spaces and is padded on the right
        info_frame.insert(END,'{:8}{:8}{:8}\n'.format(key,begin,end))
    #the description is added to the console, so the user can see the protein name and further information in the Swissprot entry
    info_frame.insert(END,description)

#here we set up the GUI, this uses a Tkinter grid layout to keep the widgets lined up
text_label = Label(window,text="Swissprot accession number: ")
text_label.grid(column=0, row=0,sticky=W)
protein_number = Entry(window,width=22)
protein_number.grid(column=1,row=0,sticky=W)
"""the send_button object is the button that runs the main update() function, which is the main function
responsible for string validation, getting the protein data, and updating the console
"""
send_button = Button(window,text="Find motifs",command = update)
send_button.grid(row=1,column=1,sticky=W)
#info_frame is the "console" which shows the sequence positions and protein descriptions
info_frame = Text(window,font=('Courier',15),width=40,height=20,wrap=WORD)
info_frame.grid(row=2,columnspan=2)
#mainloop initiates the tkinter window and runs the program
mainloop()

