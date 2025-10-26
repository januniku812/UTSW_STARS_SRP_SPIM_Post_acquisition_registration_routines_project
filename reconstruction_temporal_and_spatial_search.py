'''

NOTE: This code cannot be rerun without the imaging data from the the Bioinformatics Department at UTSW 
TEMPORAL AND SPATIAL SEARCH AND SORT ALGORITHM 
'''

# necessary imports for reading tif images, SSIM calculations, and possible mappings
import cv2 # needed for reading images in the cell directory
import os # needed for parsing file paths in the directory
from skimage.metrics import structural_similarity as ssim
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import pearsonr
import numpy as np
from tifffile import imread, imsave, TiffFile# reading in tiffiles for SSIM + stacking and then saving them later for 3D timepoint generation


# cosine generation function for later
def cosine_func(x, A, B, C, D):
    return A * np.cos(B * x + C) +  D


DATA_DIR = "[]" # this needs to be set manually
RUN_NUMBER = _ # needs to be set manually

directory_of_cells = f'/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/Cell'
repeating_index_array = []
value_for_adjacent_plane_being_too_similar_in_ssim_coefficient = 1.0
run_number = 2


for i in np.arange(1, 51): # the end parameter is set manually based on 'z positions' were imaged using SPIM microscope manually 

    reference_cell_path = f'/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/Cell' + str(i) + '/1_CH00_000000.tif'
    cell_path = f'/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/Cell' + str(i)

    reference_temporal_xy_slice = imread(reference_cell_path)
    if reference_temporal_xy_slice.shape != (2048, 2048):
        reference_temporal_xy_slice = reference_temporal_xy_slice[0]


    # saving a list of all the tif files within the cell we want to traverse
    image_paths = []

    for dirpath, dirname, filenames in os.walk(cell_path):
        for filename in filenames:
            if filename.lower().endswith(".tif"):
                image_paths.append(os.path.join(dirpath, filename))


    print(f"image paths: {image_paths}")

    pearson_values = [] # where the ssim value at index 0 should be 1 as it refers to the SSIM index between the reference image and itself as the first image in the directory
    for temp_index in np.arange(0, 200):
        if 0 <= temp_index < 10:
            temporally_spaced_xy_slice_path = cell_path + '/1_CH00_00000' + str(temp_index) + '.tif'
        elif 10 <= temp_index < 100:
            temporally_spaced_xy_slice_path = cell_path + '/1_CH00_0000' + str(temp_index) + '.tif'
        else:
            temporally_spaced_xy_slice_path = cell_path + '/1_CH00_000' + str(temp_index) + '.tif'

        temporally_spaced_xy_slice = imread(temporally_spaced_xy_slice_path)

        if temporally_spaced_xy_slice.shape != (2048, 2048):
            temporally_spaced_xy_slice = temporally_spaced_xy_slice[temp_index]


        print(f"temporally spaced xy slice path: {temporally_spaced_xy_slice_path}")
        print(f"reference xy slice path!!: {reference_cell_path}")
        fig, axs = plt.subplots(1, 2, figsize=(8,4))
        axs[0].imshow(temporally_spaced_xy_slice, cmap="gray")
        axs[1].imshow(reference_temporal_xy_slice, cmap="gray")
        plt.title(f"reference to temporal slice @ {temporally_spaced_xy_slice_path}: z index {i} temporal index {temp_index} ")
        plt.show()


        print(f"reference temporal slice: {reference_temporal_xy_slice.shape}")
        print(f"temporally spaced xy slice: {temporally_spaced_xy_slice.shape}")

        # creating a mask so that the ssim only compares the non black pixels
        reference_temporal_slice_mask = np.any(reference_temporal_xy_slice != 0, axis=-1)
        temporally_spaced_slice_mask = np.any(temporally_spaced_xy_slice != 0, axis=-1)
        reference_temporal_xy_slice = reference_temporal_xy_slice[reference_temporal_slice_mask]
        temporally_spaced_xy_slice = temporally_spaced_xy_slice[temporally_spaced_slice_mask]
       # temporally_spaced_xy_slice = cv2.resize(temporally_spaced_xy_slice,
                                             # (1024,1600) # width, height which then get flipped to ensure that the final dim
                                             # )

        print(f"reference temporal slice: {reference_temporal_xy_slice.shape}")
        print(f"temporally spaced xy slice: {temporally_spaced_xy_slice.shape}")

        pearsonr_score = round((pearsonr(reference_temporal_xy_slice.flatten(),
                          temporally_spaced_xy_slice.flatten()) # effectively avoiding a ValueError for win_size exceeds image extent (the images are technically 3 channels)
                            [0]), # ensuring that the first value of the ssim function return tuple is taken as the final ssim score
                          6)

        print(f"pearson score: {pearsonr_score}")

        pearson_values.append(pearsonr_score)

    print("Pearson Correlation Coefficient Values: " + str(pearson_values))

    largest_value_indices = [index for index, value in enumerate(pearson_values) if value == pearson_values[0]] # here the "largest" value is the pearson correlation coefficient of the first slice to itself as it marks the repeat of the cardiac cycle
    print(f"largest value indices: {largest_value_indices}")
    largest_value_index = (np.sum(largest_value_indices)/21)

    repeating_index_array.append(largest_value_index)

    repeating_cycle_info = "The slice at which the cycle repeats for cell " + str(i) + " again is : " + str(largest_value_index)
    print(repeating_cycle_info)
    with open(f"{DATA_DIR}_{RUN_NUMBER}_info.txt", "a") as f:
        f.write(repeating_cycle_info)

    # plotting the pearson correlation coefficient values
    plt.figure()
    x = np.arange(len(pearson_values)) # creating x values from 0 to the size of the data so they can form coordinates with the wave data
    y_data = np.array(pearson_values)
    plt.scatter(x, y_data, color='red', label="scatter plot of data")
    plt.plot(x, y_data, color='blue', label="continous plot of data")

    print(f"Amplitude: {(max(pearson_values) - min(pearson_values))/2}")

    # optimizing a fitted cosine wave by starting with noisy that and then fitting it
    p0 = [((max(pearson_values) - min(pearson_values)) /2), # amplitude = (max y - min y)/2
          (2 * np.pi * (6.25)) / len(x), # the b value here adheres to th enumber of data point but there should be around 250/57.1 ~ cycles
          0, # there is no initial phase shift as the first sssim value is close to 1
          np.mean(pearson_values)] # intial guess for A,B,C,D in Acos(Bx+C) + D
    noisy_func = p0[0] * np.cos((p0[1]*x) + p0[2]) + p0[3]

    params, covariance = curve_fit(cosine_func, x, noisy_func, p0=p0)
    print("the phase shift for z position " + str(i+1) + " is " + str(params[2]))
    fitted_function = cosine_func(x, params[0], params[1], params[2], params[3])
    plt.xlabel('Index')
    plt.ylabel('PCC Value with regards to first reference slice')
    plt.title('PCC values through temporally spaced xy sheets for ' + str(i))
    plt.grid(True)
    plt.show()


for index, repeating_index in enumerate(repeating_index_array):
    print("the repeating index at z position " + str(index+1) + " is " + str(repeating_index))

repeating_index_across_cells = 66 # change as needed based on the code above 
print(f"average repeating index is {repeating_index_across_cells}")

total_start = time.time()

corresponding_temporal_image_stacks = []
print(f"correspondming temporal image stack shape check: {np.array(corresponding_temporal_image_stacks).shape}")
technical_starting_indexs_corresponding_to_z_position = [0]
starting_index = 1
for z_position_index in np.arange(starting_index, starting_index+51):
    print(f"working on z position {z_position_index}")
    cell_directory_path = f'/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/Cell' + str(z_position_index) + '/'
    if z_position_index == starting_index:
        start = time.time()
        technical_starting_point = technical_starting_indexs_corresponding_to_z_position[z_position_index-starting_index]
        image_stack = []
        for within_temporal_stack_index in np.arange(0, repeating_index_across_cells):
            slice_index = technical_starting_point + within_temporal_stack_index
            image_path = cell_directory_path + "1_CH00_000" + ("0" * (3-len(str(slice_index)))) + str(slice_index)+".tif"
            image = imread(image_path, key=slice(None))
            print(f"index being added to: {z_position_index-starting_index}")
            image_stack.append(np.array(image))

        end = time.time()
        info = f"it took {end - start} seconds to establish the first z position's start index"
        file_path = f"/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/non_sliding_window_approach_run_{run_number}/{DATA_DIR}_indexing_information_{run_number}.txt"
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "a") as f:
            f.write(info)

        print(info)
        corresponding_temporal_image_stacks.append(image_stack)

    else:
        total_time_for_saving_stack_start = time.time()
        size_of_cell_directory = np.sum(np.fromiter((1 for file in os.listdir(cell_directory_path) if file.lower().endswith(".tif")), dtype=int))
        best_ssim_correlation_coefficient = 0
        best_starting_index = 0
        start = time.time()
        temporal_window_stack = []
        reference_slice = corresponding_temporal_image_stacks[z_position_index -starting_index - 1][0]
        end_point = repeating_index_across_cells
        for potential_technical_starting_index in np.arange(0, end_point):
            start = time.time()
            slice_index = potential_technical_starting_index
            print(f"slice index: {slice_index}")
            image_path = cell_directory_path + "1_CH00_000" + ("0" * (3-len(str(slice_index)))) + str(slice_index) + ".tif"
            print(image_path)

            potential_starting_slice = None
            with TiffFile(image_path) as tif:
                array_version = np.array(tif)
                img = np.resize(array_version, (1600, 1024))
                sections = []
                for page in tif.pages:
                    section = page.asarray()
                    if section.shape != (1600, 1024):
                        section = section.T
                        sections.append(section)

                potential_starting_slice = sections[0]

            print(f"potential starting slice: {potential_starting_slice.shape}")
            print(f"reference_slice: {reference_slice.shape}")
            ssim_score = ssim(np.array(potential_starting_slice), np.array(reference_slice))
            end = time.time()

            test_info_time = f"it took {end-start} seconds to test technical temporal starting index {potential_technical_starting_index} for z position: {z_position_index}"
            print(test_info_time)

            file_path = f"/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}_reconstructed_temporal_stacks/{DATA_DIR}_indexing_information_{run_number}.txt"
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            with open(file_path, "a") as f:
                f.write(test_info_time)

            if ssim_score > best_ssim_correlation_coefficient:
                best_starting_index = potential_technical_starting_index
                best_ssim_correlation_coefficient = ssim_score

        for within_temporal_window_index in (best_starting_index, best_starting_index + repeating_index_across_cells):
            within_temporal_window_index_path = cell_directory_path + "1_CH00_000" + ("0" * (3-len(str(within_temporal_window_index)))) + str(within_temporal_window_index)+".tif"
            within_temporal_window_image = imread(within_temporal_window_index_path)
            if within_temporal_window_image.shape != (2048, 2048):
                within_temporal_window_image = within_temporal_window_image[0]
            temporal_window_stack.append(within_temporal_window_image)

        corresponding_temporal_image_stacks.append(temporal_window_stack)
        technical_starting_indexs_corresponding_to_z_position.append(best_starting_index)

        total_end_for_saving_temporal_stack = time.time()
        index_finding_info = f"it took {total_end_for_saving_temporal_stack - total_time_for_saving_stack_start} seconds to find the potential technical starting index {best_starting_index} for z pos {z_position_index} and save the temporal stack for z position {z_position_index}"
        print(index_finding_info)

        file_path = f"/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}/reconstructed_temporal_stacks/{DATA_DIR}_indexing_information_{run_number}.txt"
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "a") as f:
            f.write(index_finding_info)

        imsave(f"/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}/reconstructed_temporal_stacks/temporal_stack_z_pos_{z_position_index}.tif", temporal_window_stack, bigtiff=True)

all_index_info = f"technical starting indexes: {technical_starting_indexs_corresponding_to_z_position}"
print(all_index_info)

file_path = f"/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}/indexing_information/{DATA_DIR}_indexing_information_{run_number}.txt"
os.makedirs(os.path.dirname(file_path), exist_ok=True)
with open(file_path, "a") as f:
    f.write(all_index_info)

timepoints_storage = []

print("SAVING THE TIMEPOINTS")
for timepoint_idx in np.arange(repeating_index_across_cells):
    timepoint = []
    for z_pos in np.arange(len(corresponding_temporal_image_stacks)):
        image = np.stack(corresponding_temporal_image_stacks[z_pos][timepoint_idx], axis=0)
        timepoint.append(image)

    output_path = f'/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}/non_sliding_window_approach_reconstructed_timepoints/timepoint_{timepoint_idx+1}.tif'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    imsave(output_path, timepoint.astype(np.uint16), bigtiff=True)
    timepoints_storage.append(timepoint)

total_end = time.time()
file_path = f"/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}_reconstructed_temporal_stacks/{DATA_DIR}_indexing_information_{run_number}.txt"
os.makedirs(os.path.dirname(file_path), exist_ok=True)
with open(file_path, "a") as f:
    f.write(f"the total run time for the reconstruction algorithm was {total_end - total_start} seconds")

print("SAVING REPEATING INDEX ARRAY")

file_path = f'/archive/bioinformatics/Danuser_lab/zebrafish/raw/Bo-Jui/KRDL_GFP_vas/{DATA_DIR}/sliding_window_approach_run_{run_number}/best_indexes_array/repeating_index_array_{run_number}.txt'
os.makedirs(os.path.dirname(file_path), exist_ok=True)
with open(file_path, 'a') as f:
    f.write(f'repeating index array: {repeating_index_array}')
	
