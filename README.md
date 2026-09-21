[sm2524@gadi-login-05 ~]$ mv results logs R/fa_dmr_wb/
mv: target 'R/fa_dmr_wb/' is not a directory
[sm2524@gadi-login-05 ~]$ pwd
/home/272/sm2524
[sm2524@gadi-login-05 ~]$ ls -la
total 80
drwx------   8 sm2524 vk72      4096 Sep 21 22:06 .
drwxr-xr-x 477 root   root     36864 Sep 21 17:40 ..
-rw-------   1 sm2524 vk72       393 Sep 16 22:04 .bash_history
-rw-r--r--   1 sm2524 nci-i272   141 Aug 20 21:10 .bash_profile
-rw-r--r--   1 sm2524 nci-i272  1551 Aug 20 21:10 .bashrc
drwxr-xr-x   4 sm2524 nci-i272  4096 Aug 20 21:39 .config
drwx------   2 sm2524 vk72      4096 Aug 20 21:49 .ssh
drwxr-xr-x   4 sm2524 vk72      4096 Sep 21 21:54 R
drwxr-xr-x   2 sm2524 vk72      4096 Sep 21 22:06 logs
drwxr-xr-x   3 sm2524 nci-i272  4096 Aug 20 21:39 ondemand
drwxr-xr-x   2 sm2524 vk72      4096 Sep 21 22:06 results
[sm2524@gadi-login-05 ~]$ find ~ -iname "fa_dmr_wb" -type d
[sm2524@gadi-login-05 ~]$
