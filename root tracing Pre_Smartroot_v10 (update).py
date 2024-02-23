#!/Users/Jasper/anaconda2/bin/python
import numpy as np
import cv2, os, math
from scipy import ndimage
from operator import itemgetter
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement
import pandas as pd
from itertools import islice
from scipy.spatial.distance import cdist

def Start():
    # Script to find all tiff files in directory of .py file
    global file_path
    file_path = os.path.dirname(os.path.abspath(__file__))
    File_List = PathItems(file_path)
    File_List = sorted(File_List) #List of all tif files

    #Iterating over all the tif files in folder
    for file_name in File_List:
        #Creating filenames and result folder
        print file_name

        total_list = [] #This is the list containing all root information of one .tif used to write .rsml
        image_read = cv2.imread(file_name) #loading the .tif into matrix

        make_grayscale(image_read, file_name)
        list_index, contours, Blue_channel_threshold = Contour_finding(image_read)

        Blue_channel_threshold = np.stack((Blue_channel_threshold,) * 3, axis=-1)

        #Here the program iterates over every selected contour. Final goal is to extend the total list, which is used for the
        #creation of rsml files
        for item in list_index:
            mask = np.zeros_like(image_read)  # Create mask where white is what we want, black otherwise
            cv2.drawContours(mask, contours, item, 255, -1)  # Draw filled contour in mask

            mask1 = np.zeros_like(Blue_channel_threshold)
            mask = (mask[:, :, 0] == 255)
            mask1[mask] = Blue_channel_threshold[mask]

            # Again setting up extremes, used to crop the image and analyse one root at a time
            extLeft = tuple(contours[item][contours[item][:, :, 0].argmin()][0])
            extRight = tuple(contours[item][contours[item][:, :, 0].argmax()][0])
            extTop = tuple(contours[item][contours[item][:, :, 1].argmin()][0])
            extBot = tuple(contours[item][contours[item][:, :, 1].argmax()][0])
            crop_img = mask1[extTop[1] - 2:extBot[1] + 2, extLeft[0] - 2:extRight[0] + 2]

            global Xrange_cropped, Yrange_cropped
            Xrange_cropped = len(crop_img[0])
            Yrange_cropped = len(crop_img)

            #name to save the images of every root (might be removed at a later stage)
            distance, vector_list, factor = Euclidean_distance_mapping(crop_img)

            if [extBot[1] - extTop[1] + 2, extBot[0] - extLeft[0] + 2] not in vector_list:
                vector_list.append([extBot[1] - extTop[1] + 2, extBot[0] - extLeft[0] + 2])

            #By making use of Dijkstra's algorithm (used in GPS systems) the 'shortest' path is found
            #In this case shortest is the one with the highest values in the distance map
            #Hence the extreme low values of non-root pixels. I had to run this twice, from right to left and vice versa.
            #For every pixel in the image a route is created, starting from the Top pixel value of the segmented root
            #The output is a dictonary with the best route. For every extreme of the LR and MR, the route is placed in the list
            # 'new_vector_list'.
            if len(vector_list) > 0:
                new_vector_list = []

                df = Dataframe_manager(distance)

                for i in range(len(vector_list)):
                    route = df[(df.Row == vector_list[i][0]) & (df.Column == vector_list[i][1])].Route.values
                    if len(route) > 0:
                        route = route[0]
                        if route != '0':
                            route = route[1:].strip().split()
                            vector = [[int(item.split('c')[0][1:]), int(item.split('c')[1])] for item in route]
                            vector.append(vector_list[i])
                            new_vector_list.append(vector)


                if len(new_vector_list) > 0:
                    main_root, laterals = Root_appointer(new_vector_list)
                    main_root = Selection(distance, [main_root],factor)

                    if len(main_root) > 0:
                        main_root = main_root[0]
                        MR_DL = [distance[item[0]][item[1]][0] for item in main_root]
                        MR_DL_min = min(MR_DL[int(len(MR_DL)*0.3):])

                        main_root = filter(lambda a: distance[a[0]][a[1]][0] != 10 ** 10, main_root)  # remove all empty values

                        # Two checks here, one to remove highly similar laterals (if two extremes are found in one lateral)
                        # Two routes may be created. The other to remove laterals crossing large parts of negative values
                        if len(laterals) > 0:
                            laterals = Selection(distance, laterals, factor)
                        if len(laterals) > 0:
                            laterals = lateral_check(laterals)

                        laterals = [item for item in laterals if min([distance[subitem[0]][subitem[1]][0] for subitem in item]) > MR_DL_min]
                        laterals = [filter(lambda a: distance[a[0]][a[1]][0] != 10 ** 10, item) for item in laterals]

                        total_list.append([main_root, laterals, extLeft, extTop, factor, distance])

        #Final creation of the .rsml
        XML_writer(total_list, file_name[:-4])

def PathItems(file_path): #Script to find all tiff files in directory of .py file
    return [item for item in os.listdir(file_path) if item[-3:] == 'tif']

def make_grayscale(image_read, file_name):
    SmartRootFolder = file_path + '/Improved_RSA/'
    name = SmartRootFolder + file_name + 'f'
    if not os.path.exists(SmartRootFolder):
        os.makedirs(SmartRootFolder)
    RSA_improve = image_read.astype(np.uint8)
    RSA_improve = cv2.cvtColor(RSA_improve, cv2.COLOR_BGR2GRAY)
    RSA_improve = np.invert(RSA_improve)
    RSA_improve = cv2.blur(RSA_improve, (3, 3))
    cv2.imwrite(name, RSA_improve)

def Contour_finding(image_read):
    # Extracting the blue channel
    Blue_channel = image_read
    Blue_channel[:, :, 1] = 0
    Blue_channel[:, :, 2] = 0

    # Setting an adaptive threshold. The hypocotyl is below the threshold, so image only contain roots
    # Next contours/objects are detected using the cv2 module.
    # https://docs.opencv.org/3.4.0/d9/d8b/tutorial_py_contours_hierarchy.html
    # Hierachy describes the relation between the contours

    imgray = cv2.cvtColor(Blue_channel, cv2.COLOR_BGR2GRAY)
    block_size = len(imgray) / 3
    if block_size % 2 == 0:
        block_size += 1
    Blue_channel_threshold = cv2.adaptiveThreshold(imgray, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, block_size, -2)
    none, contours, hierarchy = cv2.findContours(Blue_channel_threshold, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)

    # Specific contours are selected based on size and location, to remove noise and contours of the square plate, respectively.
    minimal_area = float(len(image_read) * len(image_read[0])) * 0.00008  # 0.000115
    maximal_area = float(len(image_read) * len(image_read[0])) * 0.3

    pixel_threshold = float(len(image_read[0])) * 0.06

    contours = [item for item in contours if (cv2.contourArea(item) > minimal_area) and (cv2.contourArea(item) < maximal_area)]

    list_index = []
    for i in range(len(contours)):
        extLeft = tuple(contours[i][contours[i][:, :, 0].argmin()][0])
        extRight = tuple(contours[i][contours[i][:, :, 0].argmax()][0])
        extTop = tuple(contours[i][contours[i][:, :, 1].argmin()][0])
        extBot = tuple(contours[i][contours[i][:, :, 1].argmax()][0])
        if extTop[1] > pixel_threshold and extTop[1] < (len(image_read) / 2):
            if extBot[1] < (len(image_read) - pixel_threshold) and extBot[1] > pixel_threshold:
                if extRight[0] > pixel_threshold * 1.5 and extLeft[0] < (len(image_read[1]) - pixel_threshold * 1.5):
                    if float((extRight[0] - extLeft[0])) / float((extBot[1] - extTop[1])) < 1:
                        list_index.append(i)

    # Selected contours are selected again based on their shape and hierarchy.
    list_index_new = []
    for item in list_index:
        # Parent of contour must not be found in contour list. This removes inner contours of touching laterals.
        # Furthermore the circumference of the minimum Enclosing Circle, must be large. This to remove selected circular
        # contours (imaging artifacts, yeast contamination and so on). Roots have huge minimum Enclosing Circles due to
        # their length
        if hierarchy[0][item][3] not in list_index:
            (x, y), radius = cv2.minEnclosingCircle(contours[item])
            if (radius ** 2) * math.pi > cv2.contourArea(contours[item]) * 5:
                list_index_new.append(item)
    list_index = list_index_new  # Final list of indexes of selected contours, normally 4 in our setup

    # (works 100% on my dataset, no detection of other stuff)
    return list_index, contours, Blue_channel_threshold

def Euclidean_distance_mapping(image_read): #Core of the program.
    image_read[image_read < 140] = 0
    image_read[image_read >= 140] = 255
    image_read = cv2.blur(image_read, (3, 3))
    image_read[image_read < 110] = 0
    image_read[:, :, 1] = image_read[:, :, 0]
    image_read[:, :, 2] = image_read[:, :, 0]

    distance_map_normal = ndimage.distance_transform_edt(image_read)

    factor = np.amax(distance_map_normal)+1

    #Increase distances of 0 to 10^10, for optimal separation using Dijkstra's algorithm (Vector_Stretching)
    distance_map_normal[:, :] = (factor - distance_map_normal[:, :])**2
    distance_map_normal[distance_map_normal == (factor - 0)**2] = 10 ** 10

    # Finding extremes within the root (tips of LR and MR). Requires optimization
    (y, x, z) = np.where(distance_map_normal == (factor - 1)**2)

    vector_init_list = []
    for i in range(0, len(x), 3):
        neighbour_list = neighbours(y[i], x[i], 'no')
        if len(filter(lambda a: distance_map_normal[a[0]][a[1]][0] == 10 ** 10, neighbour_list)) >= 5:
            vector_init_list.append([y[i], x[i]])

    return distance_map_normal,vector_init_list, factor

def routeplanner(df, iterator):
    #Loop over dataframe indexes. Iterator = all indexes to check
    for i, row in islice(iterator, 0, None):
        if i != 0:
            #Calculate all distances from POI and find neighbours
            df['distance2'] = np.sqrt(np.square(np.column_stack(((df['Row'] - df.loc[i, 'Row']).values, (df['Column'] - df.loc[i, 'Column']).values))).sum(axis=1))
            mask = df.distance2 < 1.42
            mask[i] = False
            neighbors_list = df.loc[mask]

            #Assign new values for POI. Cost = total 'euclidean cost' from route + 'euclidean cost' of current point
            try:
                NN = neighbors_list['Cost'].idxmin()
                df.loc[i, 'Cost'] = df.loc[i, 'Distance'] + df.loc[NN, 'Cost']
                df.loc[i, 'Route'] = '{} r{}c{}'.format(df.loc[NN, 'Route'], df.loc[NN, 'Row'], df.loc[NN, 'Column'])
            except:
                ''
    return df

#First run of Dijkstra's algorithm (from left to right)
def Dataframe_manager(matrix):
    #Create input for dataframe
    distance = matrix[:, :, 0].flatten()
    columns = range(len(matrix[0])) * len(matrix)
    rows = np.array([[i] * len(matrix[0]) for i in range(len(matrix))]).flatten()
    costs = [float('Inf')] * (len(matrix[0]) * len(matrix))

    #Create dataframe
    data = {'Row': rows, 'Column': columns, 'Cost': costs, 'Route': '', 'Distance': distance}
    order = ['Row', 'Column', 'Route', 'Cost', 'Distance']
    df = pd.DataFrame(data, columns=order)

    #Remove empty values
    df = df[df.Distance != 10 ** 10].reset_index()

    #Set start values
    df.loc[0, 'Cost'] = 0
    df.loc[0, 'Route'] = str(0)

    #Start determination of route
    df = routeplanner(df, df.iterrows())

    #Check to fix all non-assigned values
    mask = df.Route == ''
    extra_check = df.loc[mask]

    previous = len(extra_check)
    continuation = True
    count = 0
    while len(extra_check) != 0 and continuation == True:
        if count % 2 == 0:
            iterator = extra_check[::-1].iterrows()
        else:
            iterator = extra_check.iterrows()

        df = routeplanner(df, iterator)

        mask = df.Route == ''
        extra_check = df.loc[mask]

        if len(extra_check) == previous:
            continuation = False
        else:
            previous = len(extra_check)

        count += 1

    df = df.sort_values(['Row', 'Column'], ascending=[True, False])
    df = routeplanner(df, df.iterrows())

    df = df.sort_values(['Row', 'Column'], ascending=[True, True])

    return df

#Pixel list to analyse (First run of Dijkstra's algorithm)
def area(y,x, area):
    area_list = []
    x_max = x+area+1
    area_list = [[y-area,item] for item in range(x-area,x_max) if -1 < y - area < Yrange_cropped and -1 < item < Xrange_cropped]

    for y1 in range(y-area+1,y+area):
        if -1 < y1 < Yrange_cropped:
            if -1 < x - area < Xrange_cropped:
                area_list.append([y1, x-area])
            if -1 < x + area < Xrange_cropped:
                area_list.append([y1, x+area])

    area_list += [[y+area,item] for item in range(x-area,x_max) if -1 < y+area < Yrange_cropped and -1 < item < Xrange_cropped]
    return area_list

#Neighbours to check (First and second run of Dijkstra's algorithm)
def neighbours(y,x, include_0):
    if include_0 == 'yes':
        neighbours_list = [[y2, x2] for y2 in range(y - 1, y + 2) for x2 in range(x - 1, x + 2) if (-1 < x2 < Xrange_cropped and -1 < y2 < Yrange_cropped and (0 <= x2 < Xrange_cropped) and (0 <= y2 < Yrange_cropped))]
    else:
        neighbours_list = [[y2, x2] for y2 in range(y - 1, y + 2) for x2 in range(x - 1, x + 2) if (-1 < x2 < Xrange_cropped and -1 < y2 < Yrange_cropped and (x != x2 or y != y2) and (0 <= x2 < Xrange_cropped) and (0 <= y2 < Yrange_cropped))]
    return neighbours_list

def Root_appointer(new_vector_list):
    main_root = sorted(new_vector_list, key=len)[-1] # longest vector = MR
    new_vector_list.remove(main_root)  # Remove MR from list. Lateral remain

    main_root = np.array(main_root)
    new_vector_list = np.array(new_vector_list)

    laterals = []
    for vector in new_vector_list:
        distance = cdist(vector, main_root).min(axis=1)[::-1]
        if np.max(distance) > 4:
            lateral = vector[-np.argmin(distance):]
            if len(lateral) > 10:
                laterals.append(vector[-np.argmin(distance):])
    return main_root, laterals

#Simply counts all negative values in euclidean distance map. if less than 1% of pixels is negative, the vector is passed
def Selection(matrix, vector_list, factor):
    selection_count = [item for item in vector_list if len(filter(lambda a: matrix[a[0]][a[1]][0] == 10 ** 10, item)) <= (len(item) * 0.01)]
    selection_length = [item for item in selection_count if (len(item) >= 40)]
    selection_max = [item for item in selection_count if min([matrix[a[0]][a[1]][0] for a in item]) < ((factor - 1.5)**2) and item not in selection_length]
    return selection_length + selection_max

#Removal of similar laterals
def lateral_check(laterals):
    laterals_new = []
    check = True
    while check == True:
        laterals_new_temp = [laterals[0]]
        area_list = [laterals[0][0]]
        x_value = laterals[0][0][1]
        y_value = laterals[0][0][0]
        for i in range(1, 3):
            area_list += area(y_value, x_value, i)

        for k in range(len(laterals)):
            if 0 != k:
                if laterals[k][0] in area_list:
                    laterals_new_temp.append(laterals[k])
                else:
                    uniques = []
                    for i in range(len(laterals[0])):
                        if laterals[0][i] not in laterals[k]:
                            uniques.append(laterals[0][i])
                    if len(uniques) < float(min(len(laterals[0]), len(laterals[k])) * 0.6):
                        laterals_new_temp.append(laterals[k])

        laterals_new.append(sorted(laterals_new_temp, key=len)[-1])
        for item in laterals_new_temp:
            laterals.remove(item)
        if len(laterals) == 0:
            check = False
    return laterals_new

#Writing the rsml file
def XML_writer(total_list, file_name):
    root = Element('rsml')
    root.set('xmlns:po', 'http://www.plantontology.org/xml-dtd/po.dtd')
    metadata = SubElement(root, 'metadata')

    version = SubElement(metadata, 'version')
    version.text = '1'

    unit = SubElement(metadata, 'unit')
    unit.text = 'inch'

    unit = SubElement(metadata, 'resolution')
    unit.text = '300.0'

    unit = SubElement(metadata, 'last-modified')
    unit.text = 'today'

    unit = SubElement(metadata, 'software')
    unit.text = 'smartroot'

    unit = SubElement(metadata, 'user')
    unit.text = 'globet'

    unit = SubElement(metadata, 'file-key')
    unit.text = 'myimage'

    x = SubElement(metadata, 'property-definitions')

    list = [['diameter', 'float', 'cm'], ['length', 'float', 'cm'], ['pixel', 'float', 'none'],
            ['angle', 'float', 'degree'], ['insertion', 'float', 'cm'], ['lauz', 'float', 'cm'],
            ['lbuz', 'float', 'cm'], ['node-orientation', 'float', 'radian']]
    for i in range(len(list)):
        entry = SubElement(x, 'property-definition')
        label = SubElement(entry, 'label')
        label.text = list[i][0]
        label = SubElement(entry, 'type')
        label.text = list[i][1]
        label = SubElement(entry, 'unit')
        label.text = list[i][2]

    image = SubElement(metadata, 'image')
    label = SubElement(image, 'captured')
    label.text = 'today'
    label = SubElement(image, 'label')
    label.text = file_name

    scene = SubElement(root, 'scene')
    plant = SubElement(scene, 'plant')

    count = 0
    for roots in sorted(total_list, key=itemgetter(2)):
        main_root, laterals, extLeft, extTop, factor,distance_matrix = roots
        ID_string = file_name+'_root_'+str(count)+'_lat_x'
        label_string = 'root_' + str(count)

        main_root_xml = SubElement(plant, 'root')
        main_root_xml.set('ID', ID_string)
        main_root_xml.set('label', label_string)
        main_root_xml.set('po:accession', 'PO:0009005')

        #properties = SubElement(main_root_xml, 'properties')
        #label = SubElement(properties, 'rulerAtOrigin')
        #label.text = '0.0'

        geometry = SubElement(main_root_xml, 'geometry')
        polyline = SubElement(geometry, 'polyline')

        for i in range(1, len(main_root)-4, 10):
            item = main_root[i]
            root2 = SubElement(polyline, 'point')
            root2.set('x', str(item[1] - 2 + extLeft[0]))
            root2.set('y', str(item[0] - 2 + extTop[1]))

        item = main_root[-1]
        root2 = SubElement(polyline, 'point')
        root2.set('x', str(item[1] - 2 + extLeft[0]))
        root2.set('y', str(item[0] - 2 + extTop[1]))

        functions = SubElement(main_root_xml, 'functions')
        function1 = SubElement(functions, 'function')
        function1.set('name', 'diameter')
        function1.set('domain', 'polyline')

        for i in range(1, len(main_root)-4, 10):
            item = main_root[i]
            diameter = SubElement(function1, 'sample')
            diameter.text = str(2*float(factor-math.sqrt(distance_matrix[item[0]][item[1]][0])))

        item = main_root[-1]
        diameter = SubElement(function1, 'sample')
        diameter.text = str(2*float(factor-math.sqrt(distance_matrix[item[0]][item[1]][0])))

        for i in range(len(laterals)):
            ID_string = file_name+'root_'+str(count)+'_lat_'+str(i)
            Label_string = 'R'+str(count)+'lat_' + str(i)

            root1 = SubElement(main_root_xml, 'root')
            root1.set('ID', ID_string)
            root1.set('label', Label_string)
            root1.set('po:accession', 'PO:0009005')

            properties = SubElement(root1, 'properties')
            #label = SubElement(properties, 'rulerAtOrigin')
            #label.text = '0.0'

            geometry = SubElement(root1, 'geometry')
            polyline = SubElement(geometry, 'polyline')
            if len(laterals[i])-1 > 9:
                for j in range(1, len(laterals[i])-4, 3):
                    item = laterals[i][j]
                    root2 = SubElement(polyline, 'point')
                    root2.set('x', str(item[1] - 2 + extLeft[0]))
                    root2.set('y', str(item[0] - 2 + extTop[1]))
            else:
                item = laterals[i][1]
                root2 = SubElement(polyline, 'point')
                root2.set('x', str(item[1] - 2 + extLeft[0]))
                root2.set('y', str(item[0] - 2 + extTop[1]))

            item = laterals[i][-1]
            root2 = SubElement(polyline, 'point')
            root2.set('x', str(item[1] - 2 + extLeft[0]))
            root2.set('y', str(item[0] - 2 + extTop[1]))

            functions = SubElement(root1, 'functions')
            function = SubElement(functions, 'function')
            function.set('name', 'diameter')
            function.set('domain', 'polyline')

            for j in range(1, len(laterals[i])-4, 3):
                item = laterals[i][j]
                diameter = SubElement(function, 'sample')
                diameter.text = str(factor - math.sqrt(distance_matrix[item[0]][item[1]][0]))

            item = laterals[i][-1]
            diameter = SubElement(function, 'sample')
            diameter.text = str(factor - math.sqrt(distance_matrix[item[0]][item[1]][0]))

            annotations = SubElement(root1, 'annotations')

        count+= 1

    prettify(root,file_name)

#Creating a structured rsml
def prettify(elem,file_name):
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)

    reparsed = reparsed.toprettyxml(indent='  ')
    f = open('Improved_RSA/'+file_name+'.rsml', 'wb')
    f.write(reparsed)
    f.close()

Start() #Start of the function
