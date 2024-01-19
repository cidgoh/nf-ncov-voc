#!/usr/bin/env python

# import necessary modules
import requests
import base64
import configparser
import sys
import argparse
import os

# define a function to create a WordPress post
def create_wordpress_post(config_file, title, content):

    # read configuration file
    config = configparser.ConfigParser()
    config.read(config_file)

    # extract necessary information from the configuration file
    password = config['DEFAULT']['password']
    api_url = config['DEFAULT']['api_url']
    category_id = config['DEFAULT']['category_id']
    username = config['DEFAULT']['username']
    password = config['DEFAULT']['password']

    # encode the username and password to create a token
    credentials = username + ":" + password
    token = base64.b64encode(credentials.encode()).decode('utf-8')
    header = {'Authorization': 'Basic ' + token}
    
    # create data to be sent in the request
    table_rows = content.strip().split('\n')
    table_html = '<table style="width:100%; border-collapse: collapse;">\n'
    for row in table_rows:
        columns = row.split(':')
        table_html += f'<tr>\n<td style="border: 1px solid black; padding: 5px; font-size: 20px;">{columns[0].strip()}</td>\n<td style="border: 1px solid black; padding: 5px; font-size: 20px;">{columns[1].strip()}</td>\n</tr>\n'
    table_html += '</table>'

    data = {
        'title': title,
        'status': 'publish',
        'categories': category_id,
        'content' : table_html
    }
    
    # send a post request to the WordPress API to create a new post
    response = requests.post(api_url, headers=header, json=data)
    # print(data)
    print(response.text)
    

# if the script is run directly, parse command line arguments and call the create_wordpress_post function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='add update to virusmvp', add_help=False)
    parser.add_argument('-t', '--title',  help = "title of the update", required = True)
    parser.add_argument('-c', '--content', help = "Content of the update", required = True)
    parser.add_argument('-r', '--config', help = "Config file")
    args = parser.parse_args()

    # check if a configuration file is provided, otherwise use the default file
    if args.config is None and not os.path.isfile("config.ini"):
        print('\n'+'Error: please provide a configue file'+'\n')    
        parser.print_help()
        sys.exit()
    if args.config: 
        config_file = args.config
    else: 
        config_file = "config.ini"
    
    # read the content of the file specified in the command line argument
    with open(args.content, 'r') as file:
        content = file.read()

    # call the create_wordpress_post function with the necessary arguments
    create_wordpress_post(config_file, args.title, content )

