
import matplotlib

from urllib import urlretrieve
import os
import nltk


import numpy as np
from bs4 import BeautifulSoup

import codecs


from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import euclidean_distances
from sklearn.preprocessing import normalize


import csv
import xlrd
import re

import wikipedia

from nltk.corpus import stopwords

def review_to_words( raw_review ):
    # Function to convert a raw review to a string of words
    # The input is a single string (a raw movie review), and 
    # the output is a single string (a preprocessed movie review)
    #
    # 1. Remove HTML
    review_text = raw_review
    #
    # 2. Remove non-letters        
    letters_only = re.sub("[^a-zA-Z]", " ", review_text) 
    #
    # 3. Convert to lower case, split into individual words
    words = letters_only.lower().split()                             
    #
    # 4. In Python, searching a set is much faster than searching
    #   a list, so convert the stop words to a set
    stops = set(stopwords.words("english"))                  
    # 
    # 5. Remove stop words
    meaningful_words = [w for w in words if not w in stops]   
    #
    # 6. Join the words back into one string separated by space, 
    # and return the result.
    return( " ".join( meaningful_words ))


content = {}
content['cisplatin'] = review_to_words(wikipedia.page("Cisplatin").content);
content['carboplatin'] = review_to_words(wikipedia.page("Carboplatin").content);
content['taxol'] = review_to_words(wikipedia.page("Paclitaxel").content);
content['paclitaxel'] = review_to_words(wikipedia.page("Paclitaxel").content);
content['fluorouracil'] = review_to_words(wikipedia.page("Fluorouracil").content);
content['5fluorouracil'] = review_to_words(wikipedia.page("Fluorouracil").content);
content['5fu'] = review_to_words(wikipedia.page("Fluorouracil").content);
content['capecitabine'] = review_to_words(wikipedia.page("Capecitabine").content);
content['temodar'] = review_to_words(wikipedia.page("Temodar").content);
content['temozolomide'] = review_to_words(wikipedia.page("Temodar").content);
content['taxotere'] = review_to_words(wikipedia.page("Taxotere").content);
content['etoposide'] = review_to_words(wikipedia.page("Etoposide").content);
content['gemcitabine'] = review_to_words(wikipedia.page("Gemcitabine").content);
content['leucovorin'] = review_to_words(wikipedia.page("Leucovorin").content);
content['epirubicin'] = review_to_words(wikipedia.page("Epirubicin").content);
content['bleomycin'] = review_to_words(wikipedia.page("Bleomycin").content);

content['external'] = review_to_words(wikipedia.page("External_beam_radiotherapy").content);
content['externalbeam'] = review_to_words(wikipedia.page("External_beam_radiotherapy").content);
content['internal'] = review_to_words(wikipedia.page("Brachytherapy").content);



csvReader = csv.reader(codecs.open("topTreatments.csv", 'rU', 'utf-8'))

clean_treatment_word_strings = []
for row in csvReader:
	broken = row[1].split(";")
	if(broken[0] == 'x'):
		continue
	vals = []
	cur_str = ""
	for s in broken:
		bipart = s.split(', ')
		vals.append(bipart[1])
		cur_str += bipart[0] + " "

	for val in vals:
		cur_str += content.get(val, '')
	clean_treatment_word_strings.append(cur_str)










# Initialize the "CountVectorizer" object, which is scikit-learn's
# bag of words tool.  
vectorizer = CountVectorizer(analyzer = "word",
                             tokenizer = None,
                             preprocessor = None,
                             stop_words = None,
                             max_features = 400)

train_data_features = vectorizer.fit_transform(clean_treatment_word_strings)

train_data_features = train_data_features.toarray()
train_data_features = normalize(train_data_features, axis=1, norm='l1')

print train_data_features

similarities = euclidean_distances(train_data_features)
print similarities

np.savetxt("bagword_features.csv", train_data_features, delimiter = ',')
np.savetxt("distances.csv", similarities, delimiter = ',')

#vocab = vectorizer.get_feature_names()
#print "Vocab:"
#print vocab


