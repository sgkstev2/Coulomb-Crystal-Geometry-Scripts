# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 13:30:32 2025

@author: Kane
"""

#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from scipy.ndimage import binary_erosion
from scipy.optimize import least_squares

def find_and_rename_images(base_folder):  
    target_folder = os.path.join(base_folder, "To Run")
    os.makedirs(target_folder, exist_ok=True)

    config_pattern = re.compile(r"config(\d+)")
    data_frame = pd.DataFrame()

    for item in os.listdir(base_folder):
        item_path = os.path.join(base_folder, item)
        if not os.path.isdir(item_path) or not config_pattern.match(item):
            continue  

        formatted_number = re.findall(r"\d+", item)[0] if re.findall(r"\d+", item) else "N/A"

        # --------------------------------------Functions--------------------------------------------------------
        def extract_perimeter_ions(image_path, threshold=10):
            im = Image.open(image_path).convert('L').rotate(90, expand=True)
            arr = np.array(im)
            mask = arr > threshold
            inner = binary_erosion(mask)
            perimeter = mask ^ inner  
            y_perimeter, x_perimeter = np.where(perimeter)
            return x_perimeter, y_perimeter

        def Fit(x_data , y_data):
            def F(params, x, y):
                a,b,u,w = params
                return ((x-u)/a)**2 + ((y-w)/b)**2 - 1

            def residuals(params, x, y):
                return F(params, x, y)

            def initial_parameters(x , y):
                u = (np.max(x) + np.min(x)) / 2
                w = (np.max(y) + np.min(y)) / 2
                a = (np.max(x) - np.min(x)) / 2
                b = (np.max(y) - np.min(y)) / 2

                init = np.array([a, b, u, w])
                
                lower = np.array([0.5*a , 0.5*b , 0.9*u , 0.9*w])
                upper = np.array([1.5*a , 1.5*b , 1.1*u , 1.1*w])
                return init, lower, upper

            init, lower, upper = initial_parameters(x_data , y_data)
            result = least_squares(residuals, init, args=(x_data, y_data),
                                   bounds=(lower, upper), max_nfev=5000)
            return result.x

        def Length_Width(a, b):
            L = 2*a
            W = 2*b
            return L, W

        # ----------------------------------------------------------------------------------------------
        vrf = vend = eta = np.nan
        ion_records = []

        for file in os.listdir(item_path):
            if file == "Calcium_image.png":
                img_path = os.path.join(item_path, file)
                
                x, y = extract_perimeter_ions(img_path , threshold = 5)

                # filtering data to remove potential outliers
                z_scores_x = (x - np.mean(x)) / np.std(x)
                z_scores_y = (y - np.mean(y)) / np.std(y)
                th = 1.65
                x_low, x_high = 10, 640
                y_low, y_high = 10, 450
                
                mask = (
                        (np.abs(z_scores_x) < th) &
                        (np.abs(z_scores_y) < th) &
                        (x >= x_low) & (x <= x_high) &
                        (y >= y_low) & (y <= y_high)
                    )
                   
                x_f, y_f = x[mask], y[mask]
                
                
                
                

                # performing fit
                a, b, u, w = Fit(x_f, y_f)
                #print(a, b, u, w)
                
                # calculating length and width
                L, W = Length_Width(a , b)
                
                
            if file == 'trap.info':
                trap_file_path = os.path.join(item_path, file)
                with open(trap_file_path, "r") as f:
                    lines = [line.strip() for line in f.readlines()]

                inside_trap = False
                inside_ion = False
                trap_data = {}
                ion_records = []

                for line in lines:
                    if line.startswith('trap {'):
                        inside_trap = True
                        continue
                    if line.startswith('ionnumbers {'):
                        inside_ion = True
                        continue
                    if line == "}":
                        inside_trap = inside_ion = False
                        continue

                    if inside_trap:
                        parts = line.split()
                        if len(parts) >= 2:
                            key = parts[0].lower()
                            try:
                                value = float(parts[1])
                            except ValueError:
                                value = parts[1]
                            trap_data[key] = value

                    if inside_ion:
                        parts = line.split()
                        if len(parts) >= 2:
                            ion_name = parts[0]
                            try:
                                ion_number = float(parts[1])
                            except ValueError:
                                ion_number = parts[1]
                            ion_records.append({'Ion': ion_name, 'Number': ion_number})

                vrf = trap_data.get('vrf')
                vend = trap_data.get('vend')
                eta = trap_data.get('eta')    
                
            
                filled_row = {
                    'config #': formatted_number,
                    'Length': L,
                    'Width': W,
                    'eta': eta,
                    'Vrf': vrf,
                    'Vend': vend
                }
    
                for ion in ion_records:
                    ion_name = ion['Ion']
                    ion_number = ion['Number']
                    filled_row[ion_name] = ion_number
    
                filled_frame = pd.DataFrame([filled_row])
                data_frame = pd.concat([data_frame, filled_frame], ignore_index=True)
                    
                # plotting fit with crystal
                Y = np.linspace(w - b, w + b, 400)
                inside = 1 - ((Y - w)/b)**2
                inside[inside < 0] = 0
                sqrt_term = np.sqrt(inside)
                x1 = u +  a*sqrt_term 
                x2 = u - a*sqrt_term 
    
                plt.imshow(np.array(Image.open(img_path).convert('L').rotate(90, expand=True)), cmap='gray')
                plt.plot(x1, Y, 'r', linewidth=2, label='Fit')
                plt.plot(x2, Y, 'r', linewidth=2)
                plt.hlines(w + b, x2[-1], x1[-1], colors='r', linestyles='-')
                plt.hlines(w - b, x2[0], x1[0], colors='r', linestyles='-')
                plt.scatter(x_f, y_f, s=3, color='b', label='Data')
                plt.axis('equal')
                plt.legend(fontsize=10)
                plt.title(f"Config {formatted_number} — Fit & Data")
                plt.show()
                    
    data_frame.to_csv('compiled_vals.csv', index=False)
    
if __name__ == "__main__":
    base_directory = os.getcwd()
    find_and_rename_images(base_directory)

