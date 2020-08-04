import pandas as pd

df = pd.read_csv('res.csv')
x = df['runtimes'].std()
y = df['runtimes'].median()
print ('Average: '+ str(y))
print ('STD: '+ str(x))
