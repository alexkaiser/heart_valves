ffmpeg -t 2.43 -i heart_sounds.wav -af "afade=t=in:st=0:d=0.3,afade=t=out:st=2.3:d=0.1" heart_sounds_filtered.wav

# ffmpeg -i heart_sounds.wav -af "afade=t=in:st=0:d=0.3" -af "afade=t=out:st=2.3:d=0.1" heart_sounds_filtered.wav


# ffmpeg -t 2.3 -i heart_sounds.wav -vf agate -vf afade=t=in:st=0:d=0.1:curve:cub -vf afade=t=out:st=2.2:d=0.2 heart_sounds_filtered.wav
