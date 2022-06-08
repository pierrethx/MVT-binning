import numpy as np




for size in range(1,8):
    listofpoints=[]
    for i in range(size):
        for j in range(1):
            listofpoints.append((i,j))

    area=len(listofpoints)
    centerx = np.mean([t[0] for t in listofpoints])
    centery = np.mean([t[1] for t in listofpoints])
    radii = [np.sqrt((t[0]-centerx)**2+(t[1]-centery)**2) for t in listofpoints]

    rav = sum(radii)/area
    reff = np.sqrt(area/3.14)

    print(area)
    print(rav/reff)
    print(np.sqrt(3.14)/2)