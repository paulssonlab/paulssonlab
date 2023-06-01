# coding: utf-8

# In[17]:


from gurobipy import *

# In[18]:


m = Model()


# In[19]:


v0 = m.addVar()
v1 = m.addVar()


# In[20]:


m.update()


# In[21]:


m.addConstr(v0 - v1 <= 4)


# In[22]:


m.addConstr(v0 + v1 <= 4)


# In[23]:


m.addConstr(-0.25 * v0 - v1 <= 1)


# In[24]:


m.setObjective(v1, GRB.MAXIMIZE)


# In[25]:


m.params.outputflag = 0


# In[26]:


m.optimize()


# In[27]:


import matplotlib.pyplot as plt

# In[28]:


plt.plot([0, 4], [0, 4])


# In[29]:


plt.plot([4, 0], [0, 4])


# In[30]:


plt.plot([0, 4], [1, 2])


# In[31]:


plt.plot([v0.X], [v1.X], "ro")


# In[32]:


plt.show()


# In[47]:


"""n = Model()
#s=[]
for i in range(0,1):
    for j in range(0,1):
        s[i][j] = n.addVar(0,1,obj=1,vtype=GRB.INTEGER,name="node1")

for k in range(1,2):
    for l in range(0,2):
        s[k][l] = n.addVar(0,1,obj=1,vtype=GRB.INTEGER,name="node1")

for i in range(2,3):
    for j in range(0,5):
        s[i][j] = n.addVar(0,1,obj=1,vtype=GRB.INTEGER,name="node1" )

s11 = n.addVar()
s12 = n.addVar()
s21 = n.addVar()
s22 = n.addVar()
s23 = n.addVar()
s24 = n.addVar()
s25 = n.addVar()"""


# In[49]:


n = Model()


# In[98]:


x = [(1, "a"), (1, "b"), (2, "b"), (3, "c"), (1, "c"), (2, "c"), (4, "c"), (5, "c")]


# In[99]:


n.update()
print(len(x), x[0])


# In[100]:


c = [2, 3, 5, 4, 10, 6, 4, 5]


# In[110]:


c1 = [[12], [4, 6], [2], [1, 4, 1, 1, 1]]


# In[111]:


if sum(c1[1]) > sum(c1[0]):
    print(x[0])

if sum(c1[1]) < sum(c1[0]) and c1[1][0] < sum(c1[2]) and c1[1][1] < sum(c1[3]):
    print(x[1], x[2])

if sum(c1[1]) < sum(c1[0]) and c1[1][0] < sum(c1[2]) and c1[1][1] > sum(c1[3]):
    print(x[1], x[4], x[5], x[6], x[7])

if sum(c1[1]) < sum(c1[0]) and c1[1][0] > sum(c1[2]) and c1[1][1] < sum(c1[3]):
    print(x[2], x[3])

if sum(c1[1]) < sum(c1[0]) and c1[1][0] > sum(c1[2]) and c1[1][1] > sum(c1[3]):
    print(x[3], x[4], x[5], x[6], x[7])


# In[112]:


p = Model()


# In[174]:


c = [[14], [4, 11]]
s1 = p.addVar(0, 1, obj=sum(c[0]), vtype=GRB.INTEGER, name="node1")
s2 = p.addVar(0, 1, obj=sum(c[1]), vtype=GRB.INTEGER, name="node2")
p.update()
con1 = p.addConstr(s1 + s2 >= 1)


# c[0][0]

p.optimize()


# In[175]:


print(s1.X)


# In[158]:


s2.X


# In[159]:


s1


# In[160]:


sum(c[1])


# In[200]:


c = [
    [14],
    [4],
    [11],
    [4, 11],
    [1, 2],
    [1, 5, 6, 8],
]  # mother costs [i] , daughter costs[i]
s1_arr_not_divide = [0] * int(len(c) / 2)
s2_arr_divide = [0] * int(len(c) / 2)
# if we give mothers and daughters
for i in range(0, int(len(c) / 2)):
    s1 = p.addVar(0, 1, obj=sum(c[0 + i]), vtype=GRB.INTEGER, name="node1")
    s2 = p.addVar(
        0, 1, obj=sum(c[int(len(c) / 2) + i]), vtype=GRB.INTEGER, name="node2"
    )

    p.update()
    con1 = p.addConstr(s1 + s2 >= 1)
    # con2 = p.addConstr()
    p.optimize()
    print("nodes", s1.X, s2.X)
    s1_arr_not_divide[i] = s1.X
    s2_arr_divide[i] = s2.X


# In[201]:


s1_arr_not_divide


# In[202]:


s2_arr_divide
