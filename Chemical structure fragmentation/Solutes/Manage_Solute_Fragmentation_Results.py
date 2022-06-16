#!/usr/bin/env python
# coding: utf-8

# In[ ]:


Solute_Fragmentation_Results = []
with open('structures_DB_combined_fragmentation_with_pattern_sorting_results.log') as f:
    for line in f.readlines():
        Solute_Fragmentation_Results.append(line)
print(Solute_Fragmentation_Results)


# In[ ]:


r = []
while Solute_Fragmentation_Results:
    a = Solute_Fragmentation_Results.index('Fragmentation from the algorithm:\n')
    b = Solute_Fragmentation_Results[a+1:].index('\n')+a+1
    c = Solute_Fragmentation_Results[a+1:b]
    r.append(c)
    Solute_Fragmentation_Results = Solute_Fragmentation_Results[b+3:]
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


f= open("Solute_Fragmentation_Results.txt","w+")
for i in results:
    f.write(str(i)+'\n')
f.close()

