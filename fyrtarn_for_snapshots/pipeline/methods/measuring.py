import os
import warnings
import itertools
import cupy as cp
import numpy as np

from scipy.ndimage import convolve
from .utilities import proc, read_h5, where, mkdir, to_fmt, mapping

warnings.filterwarnings("ignore", message="loadtxt: input contained no data")


@proc
def measuring_proc(nucleus, directory, c):
    """
    Process job for measuring features in nuclei.

    Parameters
    ----------
    nucleus : int
        Index for the tracked nucleus.
    directory : str
        Path to the directory for input/outputs.
    c : dict
        Dictionary for config.
    """
    # set CSV format
    fmt = ["%d"] + ["%.3f"] + ["%d"] + 12*["%.3f"] + 2*["%d"] + ["%.3f"]
    header = "frame,time (s),label,ch0_z_sub (px),ch0_y_sub (px),ch0_x_sub (px)," + \
             "ch0_z (nm),ch0_y (nm),ch0_x (nm),ch1_z_sub (px),ch1_y_sub (px)," + \
             "ch1_x_sub (px),ch1_z (nm),ch1_y (nm),ch1_x (nm),sisters,replicated,distance (nm)"
    
    # store tracks and get pixel size
    tracks = []  
    pixel = c["deskewing"]["pixel"]
    
    # get paths
    frames_path = os.path.join(directory, "frames")
    nuclei_path = os.path.join(directory, "nuclei")
    corrections_path = os.path.join(directory, "corrections")
    tracks_path = mkdir(os.path.join(nuclei_path, to_fmt(nucleus), "tracks"))

    # read aberration weights
    weights = np.loadtxt(os.path.join(corrections_path, "aberrations",
                         "000-001.csv"), delimiter=",", ndmin=2)

    # get path to nucleus and read track
    nucleus_path = os.path.join(nuclei_path, to_fmt(nucleus))
    nucleus_track = np.loadtxt(os.path.join(nucleus_path, "info.csv"),
                               delimiter=",", skiprows=1, ndmin=2)
    # read channel tracks
    chnls, chnl_tracks = [0,1], []
    chnls_path = os.path.join(nucleus_path, "channels")
    for i in range(2):
        chnl = chnls[i]
        chnl_tracks.append([])

        # find in channel directory
        chnl_tracks_path = os.path.join(chnls_path, to_fmt(chnl), "tracks")
        if os.path.exists(chnl_tracks_path):
            for item in os.listdir(chnl_tracks_path):
                if item[:-4].isnumeric():
                    track = np.loadtxt(os.path.join(chnl_tracks_path, item),
                                       delimiter=",", skiprows=1, ndmin=2)
                    chnl_tracks[i].append(track)

    # get number of tracks in each channel
    counts = [len(chnl_tracks[0]), len(chnl_tracks[1])]
    if min(counts) == 0:
        return tracks  # if none in one channel
    scores = np.zeros(counts, dtype=int)
    for ij in itertools.product(range(counts[0]), range(counts[1])):
        # make pair
        i, j = ij
        pair = [chnl_tracks[0][i], chnl_tracks[1][j]]
                
        # get start and end of pair
        start = int(min(pair[0][0,0], pair[1][0,0]))
        end = int(max(pair[0][-1,0], pair[1][-1,0]))
    
        # compute score for pair
        positions = np.empty((2, 3))  
        indices, exists = np.zeros(2, dtype=int), np.zeros(2, dtype=bool)
        for frame in range(start, end+1):
            # update positions at frame (if exist)
            for k in range(2):
                if (indices[k] < len(pair[k])) and (int(pair[k][indices[k],0]) == frame):
                    positions[k] = pair[k][indices[k],2:5]
                    indices[k] += 1
                    exists[k] = True

            # increase if frames match and in radius
            search = c["measuring"]["search"]
            if np.all(exists) and (np.linalg.norm(positions[0] - positions[1]) < search):
                scores[i,j] += 1
            exists[:] = 0

    # make groups of channel track "pairs"
    groups, place = [], [[],[]]  # hold indices
    for i in range(2):
        for j in range(counts[i]):
            group = [[],[]]
            group[i].append(chnl_tracks[i][j])
            groups.append(group)
            place[i].append(j + i*counts[0])

    while True:
        # get best score
        if np.max(scores) < c["measuring"]["overlap"]:
            break  # break on minimum overlap
        i, j = np.unravel_index(np.argmax(scores), scores.shape)
        scores[i,j] = 0
        
        # arrange indices for placing channel tracks
        after, before = place[0][i], place[1][j]
        if before == after:
            continue  # already in same group
        if after > before:
            before, after = after, before

        # merge channel tracks group
        group = groups.pop(before)
        for k in range(2):
            # update indices then place in new group
            for l in range(len(place[k])):
                if place[k][l] == before:
                    place[k][l] = after
                elif place[k][l] > before:
                    place[k][l] -= 1
                else:
                    pass
            for l in range(len(group[k])):
                groups[after][k].append(group[k][l])

    # make paired tracks by piecing groups together
    for group in groups:
        # ignore if no pairs
        if len(group[0])*len(group[1]) == 0:
            continue

        # find the start/end in each group
        start, end = None, None
        for i in range(2):
            # sort for later
            group[i].sort(key=len, reverse=True)
            for track in group[i]:
                if start is None:
                    start = int(track[0,0])
                else:
                    if track[0,0] < start:
                        start = int(track[0,0])
                if end is None:
                    end = int(track[-1,0])
                else:
                    if track[-1,0] > end:
                        end = int(track[-1,0])

        # build lines for each track
        info, exists = [], np.zeros(2, dtype=bool)
        for frame in range(start, end+1):
            # get subset lower bound for frame coordinates
            nucleus_idx = where(nucleus_track[:,0] == frame)
            label = int(nucleus_track[nucleus_idx, 1])
            lower = nucleus_track[nucleus_idx, 5:8]
            
            # initialize line
            line = np.empty(18)
            line[0], line[1] = frame, c["general"]["time"]*frame  # time
            line[2] = label
            line[3:] = np.nan
            line[15:17] = 0  # sisters/replicated

            # fill in line
            for i in range(2):
                # use longest tracks last
                for track in group[i]:
                    idx = where(track[:,0] == frame)
                    if (idx is not None) and (not np.isnan(track[idx,5:8]).any()):
                        # mark that localization exists for frame
                        exists[i] = True
                        
                        # add positions and sisters
                        line[3+i*6:6+i*6] = track[idx,5:8]
                        line[6+i*6:9+i*6] = lower + track[idx,5:8]
                        line[15] += track[idx,8]  # sisters
                        break
            
            # ch0 localization exists
            if exists[0]:
                # get position (in nm) of ch0 in ch1
                chromatic = mapping(line[6:9]) @ weights
                line[6:9] = pixel * (line[6:9] - chromatic)

            # ch1 localization exists
            if exists[1]:
                line[12:15] *= pixel  # position (in nm)

            # add distance
            if exists[0] and exists[1]:
                line[17] = np.linalg.norm(line[6:9] - line[12:15])
                
            # clear and append
            exists[:] = 0
            info.append(line)

        # get track as array and esimate replication start with window
        info = np.asarray(info)
        idx = where(convolve(info[:,15], np.ones(5), mode="nearest") >= 3)
        if idx is not None:
            info[idx:,16] = 1
        
        # save as CSV
        np.savetxt(os.path.join(tracks_path, to_fmt(len(tracks)) + ".csv"), info,
                   fmt=fmt, header=header, delimiter=",", comments="")
        tracks.append(info)  # add to tracks
    
    return tracks
