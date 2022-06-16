#!/usr/bin/env python
# coding: utf-8

# In[ ]:


IL_Fragmentation_Results = []
with open('structures_DB_combined_fragmentation_with_pattern_sorting_results.log') as f:
    for line in f.readlines():
        IL_Fragmentation_Results.append(line)


# In[ ]:


print(IL_Fragmentation_Results)


# In[ ]:


r = []
while IL_Fragmentation_Results:
    a = IL_Fragmentation_Results.index('Fragmentation from the algorithm:\n')
    b = IL_Fragmentation_Results[a+1:].index('\n')+a+1
    c = IL_Fragmentation_Results[a+1:b]
    r.append(c)
    IL_Fragmentation_Results = IL_Fragmentation_Results[b+3:]
print(r)
print(len(r))


# In[ ]:


results = []
for i in range(0,len(r)):
    results.append([0]*122)
print(len(results)) 

for i in range(0,len(r)):
    for j in r[i]:
        d = j.split()
        results[i][int(d[-2])-1] = int(d[-1])
print(results)


# In[ ]:


f= open("IL_Fragmentation_Results.txt","w+")
for i in results:
    f.write(str(i)+'\n')
f.close()

