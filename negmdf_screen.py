import os
import csv
import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import Point, Polygon, LineString

# Extract sample data
def sample_data_generator(sample_path):
    sample_data = []
    with open(sample_path, 'r', newline='', encoding='utf-8') as f:
        for row in f:
            if row.split(',')[0] == 'id':
                pass
            else:
                list_tem = [
                    row.split(',')[0],
                    row.split(',')[1],
                    row.split(',')[2],
                    int(row.split(',')[3]),
                    float(row.split(',')[4].strip())
                ]
                sample_data.append(list_tem)
    return sample_data

# Extract compound names
def compounds_name_generator(target_compound_path):
    compounds_name = []
    with open(target_compound_path, 'r', newline='', encoding='utf-8') as f:
        for row in f:
            if row.split(',')[0] == 'compound':
                pass
            else:
                compounds_name.append(row.split(',')[0].strip())
    return compounds_name

# Extract compound features
def compounds_feature_generator(target_compound_path):
    with open(target_compound_path, 'r', newline='', encoding='utf-8') as f:
        compounds_feature = {}
        for row in f:
            compound_feature = []
            if row.split(',')[0] == 'compound':
                pass
            else:
                i = 3
                while i + 1 < len(row.split(',')):
                    if row.split(',')[i].strip() == '':
                        break
                    else:
                        compound_feature.append([
                            int(row.split(',')[i].strip()),
                            float(row.split(',')[i + 1].strip())
                        ])
                    i += 2
                compounds_feature[(float(row.split(',')[2]), row.split(',')[1])] = compound_feature
    return compounds_feature

# Generate data for convex hull
def hull_data_generator(initial_mcr, compound):
    compound_arrays = []
    for array in compound:
        new_array = []
        for i in range(array[0] + 1):
            new_array.append(i * array[1])
        compound_arrays.append(new_array)

    cartesian_product = [[]]
    for array in compound_arrays:
        new_result = []
        for x in cartesian_product:
            for y in array:
                new_result.append(x + [y])
        cartesian_product = new_result
        

    hull_data = []
    for item in cartesian_product:
        mcr = initial_mcr
        for num in item:
            mcr += num
        hull_data.append([int(mcr), mcr - int(mcr)])
    
    return hull_data


# Screen for one compound
def single_screening(compounds_feature, sample_data, tolerance=2e-2):
    for key in compounds_feature:
        # Initialize screening result with headers
        sample_data_screening = [['id', 'rt', 'mz', 'Integer', 'decimal']]
        
        # Generate convex hull data
        hull_data = hull_data_generator(key[0], compounds_feature[key])
        hull_array = np.array(hull_data)
        hull = ConvexHull(hull_array)
        convex_hull_points = hull_array[hull.vertices]
        convex_hull_polygon = Polygon(convex_hull_points)
        
        # Process each sample point
        for point in sample_data:
            if int(point[3] % 2) == int(key[1]):
                test_point = Point(point[3], point[4])
                test_point_shapely = Point(test_point)
                
                # Check if point is inside or touching the polygon
                is_inside = convex_hull_polygon.contains(test_point)
                is_touching_edge = convex_hull_polygon.touches(test_point)
                
                if is_inside or is_touching_edge:
                    sample_data_screening.append(point)
                else:
                    # Check proximity to polygon edges
                    check = 0
                    convex_hull_vertices = np.array(convex_hull_points)
                    for i in range(len(convex_hull_vertices)):
                        p1 = convex_hull_vertices[i]
                        p2 = convex_hull_vertices[(i + 1) % len(convex_hull_vertices)]
                        line = LineString([p1, p2])
                        distance = test_point_shapely.distance(line)
                        if distance < tolerance:
                            check = 1
                            break
                    
                    if check == 1:
                        sample_data_screening.append(point)
        
        # Return the screened data for the compound
        return sample_data_screening

# Screen for multiple compounds
def multiple_screening(mode, sample_path, compounds_feature_path, output):
    # Generate compound names and features
    compounds_name = compounds_name_generator(compounds_feature_path)
    compounds_feature = compounds_feature_generator(compounds_feature_path)

    if mode == "single":
        screening_result = []
        sample_points = sample_data_generator(sample_path)
        i = 0
        for key in compounds_feature:
            screening_result.append([compounds_name[i]])
            screening_result += single_screening({key: compounds_feature[key]}, sample_points)
            i += 1
        
        with open(output, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            for row in screening_result:
                writer.writerow(row)


    elif mode == "multiple":
    
        # Traverse sample path
        for root, _, files in os.walk(sample_path):
            for file in files:
                file_path = os.path.join(root, file)
                screening_result = []

                # Only process files containing 'new' in their name
                if file_path.endswith(".csv"):
                    sample_points = sample_data_generator(file_path)
                    i = 0
                    
                    # Screen each compound feature
                    for key in compounds_feature:
                        screening_result.append([compounds_name[i]])
                        screening_result += single_screening({key: compounds_feature[key]}, sample_points)
                        i += 1

                    # Generate output file path
                    screening_result_path = os.path.join(
                        output, os.path.basename(file_path).replace(".csv","_screened.csv")
                    )

                    # Write screening results to CSV
                    with open(screening_result_path, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        for row in screening_result:
                            writer.writerow(row)

import argparse


def parse_arguments():
    """
    Parse command-line arguments with two subcommands: single and multiple.

    Subcommand 1: single
        - Input a CSV file containing ion list
        - Input a NegMDF window CSV file
        - Output a CSV file path

    Subcommand 2: multiple
        - Input a folder containing multiple ion list CSV files
        - Input a NegMDF window CSV file
        - Output folder path for the resulting CSV files
    """
    parser = argparse.ArgumentParser(
        description="Screen ion list based on NegMDF window."
    )
    
    # 创建子命令解析器
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommands: single or multiple")

    # 子命令1: single
    single_parser = subparsers.add_parser(
        "single", 
        help="Process a single ion list CSV file."
    )
    single_parser.add_argument(
        "-i", "--ion_list",
        type=str,
        required=True,
        help="Path to the input ion list file (CSV format)."
    )
    single_parser.add_argument(
        "-w", "--window",
        type=str,
        required=True,
        help="Path to the NegMDF window CSV file."
    )
    single_parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Path to save the output CSV file."
    )

    # 子命令2: multiple
    multiple_parser = subparsers.add_parser(
        "multiple", 
        help="Process multiple ion list CSV files in a folder."
    )
    multiple_parser.add_argument(
        "-i", "--ion_list",
        type=str,
        required=True,
        help="Path to the folder containing ion list CSV files."
    )
    multiple_parser.add_argument(
        "-w", "--window",
        type=str,
        required=True,
        help="Path to the NegMDF window CSV file."
    )
    multiple_parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Path to the folder where output CSV files will be saved."
    )

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    multiple_screening(args.command, args.ion_list, args.window, args.output)
    # if args.command == "single":
    #     multiple_screening(mode="single")
    # elif args.command == "multiple":
    #     ...
    # sample_data = sample_data_generator("data/ion_list.csv")
    # # compounds_name = compounds_name_generator("data/NegMDF_window.csv")
    # compounds_feature = compounds_feature_generator("data/NegMDF_window.csv")
    # # sample_data_screening = single_screening(compounds_feature, sample_data)
    # multiple_screening('data', 'data/NegMDF_window.csv')
    
    

    
