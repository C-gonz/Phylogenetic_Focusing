#! /usr/bin/env python3
# Script for converting space-delimited items into tab-delimited.
def delim_convert():
   try:
        with open("./hmmr_hits_table.txt", "r") as read_file:
            with open("./hmmr_hits_table.tsv", "a") as write_file:
                for line in read_file:
                    final_string = ""
                    line_list = line.split()
                    for item in line_list:
                        # print("The current item is: ", item) FOR DEBUG. LIST IS WORKING.
                        item_string = str(item)
                        item_string += "\t"
                        final_string += item_string
                        # print("Final_string is currently: ", final_string) FOR DBUG. FINAL_STRING WORKS.
                    final_string += "\n"
                    write_file.write(final_string)
   except IOError:
       print("IOError.")
delim_convert()