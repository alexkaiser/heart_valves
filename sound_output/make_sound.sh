
ffmpeg -i mitral_1202357_512_git_3735567_two_leaflet_eight_connectors_real_time.mp4 \
    -i heart_sounds_and_silence.wav -c:v libx264 -strict -2 -map 0:v -map 1:a mitral_1202357_512_git_3735567_two_leaflet_eight_connectors_real_time_SOUND.mp4


#ffmpeg -i mitral_1092466_256_git_4dc23b8_two_leaflet_updated_damping_damping_2x_real_time.mp4 \
#    -i heart_sounds_and_silence.wav -c:v libx264 -strict -2 -map 0:v -map 1:a mitral_1092466_256_git_4dc23b8_two_leaflet_updated_damping_damping_2x_real_time_SOUND.mp4


# -c:a libvorbis
# -af "volume=10.0"

# works great, some annoying clicks
# ffmpeg -i mitral_cycle_PERIODIC_9004892_git_d7914ba6_more_wrap_256_radial_real_time.mp4 \
#        -t 2.3 -i heart_sounds.wav -c:v libx264 -strict -2 -map 0:v -map 1:a mitral_cycle_PERIODIC_9004892_git_d7914ba6_more_wrap_256_radial_real_time_SOUND.mp4

