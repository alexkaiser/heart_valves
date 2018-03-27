import os 

def update_bounding_box(file_name, box = None):
	'''
	Looks for a line of the form 
	    %%BoundingBox:    xmin ymin xmax ymax
	in file_name

	Adjusts current box to take the largest available box 
	'''

	file = open(file_name,'r')
	found = False 
	for line in file:
		if line.startswith('%%BoundingBox:'): 
			found = True 
			split = line.split()

			xmin = int(split[1])
			ymin = int(split[2])
			xmax = int(split[3])
			ymax = int(split[4])
			break 

	file.close()

	assert found 

	if box is None:
		box_new = [xmin, ymin, xmax, ymax]
	else:
		xmin = min(xmin, box[0]) 
		ymin = min(ymin, box[1]) 
		xmax = max(xmax, box[2]) 
		ymax = max(ymax, box[3]) 
		box_new = [xmin, ymin, xmax, ymax]

	return box_new


def crop(file_name, box, file_name_new = None):
	'''
	Manually crop eps to specified 
	'''

	if file_name_new is None:
		if file_name.endswith('_uncropped.eps'):
			base_name = file_name.split('_uncropped.eps')[0]
			file_name_new = base_name + '.eps'
		else:
			assert False  

	file = open(file_name, 'r')
	new_file = open(file_name_new, 'w')
	
	for line in file:
		if line.startswith('%%BoundingBox:'): 
			new_file.write('%%BoundingBox:' + str(box[0]) + ' ' + str(box[1]) + ' ' + str(box[2]) + ' ' + str(box[3]) + '\n')
		else:
			new_file.write(line)

	file.close()
	new_file.close()



if __name__ == '__main__':

	one_family_plots = ['anterior_tension_plot_circ_uncropped.eps',
						'anterior_tension_plot_radial_uncropped.eps',
						'posterior_tension_plot_circ_uncropped.eps',
						'posterior_tension_plot_radial_uncropped.eps']

	box = None
	for plot in one_family_plots:
		box = update_bounding_box(plot, box)

	for plot in one_family_plots:
		crop(plot, box)






