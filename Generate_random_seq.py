#!/usr/bin/env python
# coding: utf-8

# In[5]:


from random import choice
def String(length):
    DNA=""
    for count in range(length):
        DNA +=choice("CGTA")
    return DNA


# In[28]:


for i in range(10):
    print('>Random_Seq'+str(i+1))
    print(String(250))


# In[18]:


protein_seq = 'RHKDESTNQCGPAVILMFYW'
len(protein_seq)


# In[ ]:




