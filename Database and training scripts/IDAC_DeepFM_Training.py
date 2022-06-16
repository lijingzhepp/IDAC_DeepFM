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




# 1.Import data
data = pd.read_csv('./Data_01032022_Shuffled.csv')
print(data)
sparse_features = list(data.columns)[9:]
dense_features = list(data.columns)[1:9]
target = ['column1']

print(dense_features)
print(sparse_features)




# 2.Label Encoding for sparse features,and do simple Transformation for dense features
for feat in sparse_features:
    lbe = LabelEncoder()
    data[feat] = lbe.fit_transform(data[feat])

mms = MinMaxScaler(feature_range=(0, 1))
data[dense_features] = mms.fit_transform(data[dense_features])




# 3.Count unique features for each sparse field
embedding_list = [4, 5, 6, 7]
e = embedding_list[1]

fixlen_feature_columns = [SparseFeat(feat, data[feat].max() + 1, embedding_dim=e)
                          for feat in sparse_features] + [DenseFeat(feat, 1, ) for feat in dense_features]

linear_feature_columns = fixlen_feature_columns
dnn_feature_columns = fixlen_feature_columns
feature_names = get_feature_names(linear_feature_columns + dnn_feature_columns)




# 4.Generate input data for model
train_and_validation = data.loc[0:47152,:]
test = data.loc[47153:52371,:]

validation_list = [list(i for i in range(0,5239)), 
                   list(i for i in range(5239,10478)), 
                   list(i for i in range(10478,15717)), 
                   list(i for i in range(15717,20956)), 
                   list(i for i in range(20956,26195)), 
                   list(i for i in range(26195,31434)), 
                   list(i for i in range(31434,36673)), 
                   list(i for i in range(36673,41913)), 
                   list(i for i in range(41913,47153))]
train_list = [list(i for i in range(5239,47153)), 
              list(i for i in range(0,5239))+list(i for i in range(10478,47153)), 
              list(i for i in range(0,10478))+list(i for i in range(15717,47153)), 
              list(i for i in range(0,15717))+list(i for i in range(20956,47153)), 
              list(i for i in range(0,20956))+list(i for i in range(26195,47153)), 
              list(i for i in range(0,26195))+list(i for i in range(31434,47153)), 
              list(i for i in range(0,31434))+list(i for i in range(36673,47153)), 
              list(i for i in range(0,36673))+list(i for i in range(41913,47153)), 
              list(i for i in range(0,41913))]




# 5.Create data storage

train_mse_res = []
validation_mse_res = []
test_mse_res = []

train_r2_res = []
validation_r2_res = []
test_r2_res = []

train_mae_res = []
validation_mae_res = []
test_mae_res = []

train_aard_res = []
validation_aard_res = []
test_aard_res = []




# 6.Cross-validation
dnn_structure_list = [(128, 128), (256, 256), (128, 128, 128),(64, 128, 256),
                      (256, 128, 64), (256, 256, 64), (256, 256, 256), (512,256,256),
                      (128, 128, 128, 128), (256,256,256,64), (512,512,256), (512,512,512), (512, 512, 256, 64)]

learningrate_list = [1e-04, 5e-04, 1e-03, 5e-03, 1e-02]
batchsize_list = [64, 128, 256, 512]
dropout_list = [0, 1e-02, 5e-02, 1e-01, 2e-01]
l2linear_list = [0, 1e-06, 1e-05, 1e-04]
l2embedding_list = [0, 1e-06, 1e-05, 1e-04]
l2dnn_list = [0, 1e-06, 1e-05, 1e-04]

x = dnn_structure_list[7]
y = learningrate_list[2]
z = batchsize_list[2]
a = dropout_list[0]
b = l2linear_list[2]
c = l2embedding_list[2]
d = l2dnn_list[0]

for i in range(0,9):
    for j in range (0,5):
        train = train_and_validation.loc[train_list[i],:]
        validation = train_and_validation.loc[validation_list[i],:]
        train_model_input = {name: train[name].values for name in feature_names}
        validation_model_input = {name: validation[name].values for name in feature_names}
        test_model_input = {name: test[name].values for name in feature_names}

        # 6.1 Define Model,train,predict and evaluate model
        model = DeepFM(linear_feature_columns, dnn_feature_columns, dnn_hidden_units=x, dnn_dropout=a, 
                       l2_reg_linear=b, l2_reg_embedding=c, l2_reg_dnn=d, task='regression')
        
        model.compile(Adam(y), "mse", metrics=['mse'], )
        es = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

        history = model.fit(train_model_input, train[target].values,
                            batch_size=z, epochs=100, verbose=1, 
                            validation_data=(validation_model_input, validation[target].values), callbacks=[es])

        history_name = 'DeepFM_dnn'+str(x)+'_embedding'+str(e)+'_lr'+str(y)+'_bs'+str(z)+'_drop'+str(a)+'_l2l'+str(b)+'_l2e'+str(c)+'_l2d'+str(d)+'_cv'+ str(i+1)+'_rep'+str(j+1)+'_history.npy'
        np.save(history_name, history.history)

        filename = 'DeepFM_dnn'+str(x)+'_embedding'+str(e)+'_lr'+str(y)+'_bs'+str(z)+'_drop'+str(a)+'_l2l'+str(b)+'_l2e'+str(c)+'_l2d'+str(d)+'_cv'+ str(i+1)+'_rep'+str(j+1)+'_model.h5'
        save_model(model, filename)

        train_ans = model.predict(train_model_input, batch_size=z)
        validation_ans = model.predict(validation_model_input, batch_size=z)
        test_ans = model.predict(test_model_input, batch_size=z)

        train_mse_res.append(mean_squared_error(train[target].values, train_ans))
        validation_mse_res.append(mean_squared_error(validation[target].values, validation_ans))
        test_mse_res.append(mean_squared_error(test[target].values, test_ans))

        train_r2_res.append(r2_score(train[target].values, train_ans))
        validation_r2_res.append(r2_score(validation[target].values, validation_ans))
        test_r2_res.append(r2_score(test[target].values, test_ans))

        train_mae_res.append(mean_absolute_error(train[target].values, train_ans))
        validation_mae_res.append(mean_absolute_error(validation[target].values, validation_ans))
        test_mae_res.append(mean_absolute_error(test[target].values, test_ans))

        train_aard_res.append(mean_absolute_percentage_error(exp(train[target].values), exp(train_ans)))
        validation_aard_res.append(mean_absolute_percentage_error(exp(validation[target].values), exp(validation_ans)))
        test_aard_res.append(mean_absolute_percentage_error(exp(test[target].values), exp(test_ans)))




# 7.Check the results visuaslly
print(train_mse_res)
print(validation_mse_res)
print(test_mse_res)

print(train_r2_res)
print(validation_r2_res)
print(test_r2_res)

print(train_mae_res)
print(validation_mae_res)
print(test_mae_res)

print(train_aard_res)
print(validation_aard_res)
print(test_aard_res)




# 8. Save the results in xlsx file
workbook = xlsxwriter.Workbook('IDAC_DeepFM_Continuous_DNN7_2.xlsx')
worksheet = workbook.add_worksheet()
row = 0
for i in range(0,4):
    worksheet.write (row, 0, train_mse_res[i])
    worksheet.write (row, 1, validation_mse_res[i])
    worksheet.write (row, 2, test_mse_res[i])
    worksheet.write (row, 3, train_r2_res[i])
    worksheet.write (row, 4, validation_r2_res[i])
    worksheet.write (row, 5, test_r2_res[i]) 
    worksheet.write (row, 6, train_mae_res[i])
    worksheet.write (row, 7, validation_mae_res[i])
    worksheet.write (row, 8, test_mae_res[i])
    worksheet.write (row, 9, train_aard_res[i])
    worksheet.write (row, 10, validation_aard_res[i])
    worksheet.write (row, 11, test_aard_res[i]) 
    row += 1
workbook.close()


#