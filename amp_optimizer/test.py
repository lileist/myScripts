
"""
def tryAgain(retries=0):
    if retries > 10: return
    try:
        print retries
        print succ
    except:
        retries+=1
        if retries > 5:
           succ = True
        tryAgain(retries)
""" 


retries=0
while True:
    try:
        retries +=1
        print
    except SomeSpecificException:
        continue
    break
   

#tryAgain()
