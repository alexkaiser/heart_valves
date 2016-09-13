import os 

if __name__ == '__main__':

    count = 1
    n = 2 
    while True: 
        if os.path.isfile('../mitral_tree_' + str(n) + '.vertex'): 
            print 'mitral_tree_' + str(n)    
            break  
        
        n *= 2
        count += 1
        if count > 20:
            break 
