FROM repbioinfo/testpersara
ENV DEBIAN_FRONTEND noninteractive
EXPOSE 8787
ENTRYPOINT ["tail"]
CMD ["-f","/dev/null"]
