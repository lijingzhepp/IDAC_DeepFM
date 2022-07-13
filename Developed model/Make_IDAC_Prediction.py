#!/usr/bin/env python
# coding: utf-8


# import libraries
import pandas as pd
import numpy as np
from numpy import exp
from pandas import cut
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn import metrics
from tensorflow.keras.optimizers import Adam, Adagrad
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error, mean_absolute_percentage_error
from tensorflow.keras.models import  save_model,load_model
from deepctr.layers import custom_objects
from deepctr.models import DeepFM
from deepctr.feature_column import SparseFeat, DenseFeat, get_feature_names
import xlsxwriter




# 1. Read sample data
data = pd.read_csv('./Sample.csv')
sparse_features = list(data.columns)[9:]
dense_features = list(data.columns)[1:9]
target = ['column1']




# 2. Read model
filename = 'DeepFM_dnn(512, 256, 256)_embedding5_lr0.001_bs256_drop0_l2l1e-05_l2e1e-05_l2d0_cv8_rep1_model.h5'
model = load_model(filename,custom_objects)




# 3.Label Encoding for sparse features,and do simple Transformation for dense features
for feat in sparse_features:
    lbe = LabelEncoder()
    data[feat] = lbe.fit_transform(data[feat])

mms = MinMaxScaler(feature_range=(0, 1))
data[dense_features] = mms.fit_transform(data[dense_features])




# 4.count #unique features for each sparse field
embedding_list = [4, 5, 6, 7]
e = embedding_list[1]

fixlen_feature_columns = [SparseFeat(feat, data[feat].max() + 1, embedding_dim=e)
                          for feat in sparse_features] + [DenseFeat(feat, 1, ) for feat in dense_features]

linear_feature_columns = fixlen_feature_columns
dnn_feature_columns = fixlen_feature_columns
feature_names = get_feature_names(linear_feature_columns + dnn_feature_columns)




# 5.Calculate predicted IDAC
sample = data.loc[:,:]
sample_model_input = {name: sample[name].values for name in feature_names}
sample_ans = model.predict(sample_model_input, batch_size=256)




# 6.Write predicted IDAC into excel file
workbook = xlsxwriter.Workbook('Predicted_IDAC.xlsx')
worksheet = workbook.add_worksheet()
row = 0
for i in range(0,52372):
    worksheet.write (row, 0, sample_ans[i])
    row += 1
workbook.close()






