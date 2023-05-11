# Sun Calculator for Go

[Sun Calculator for Python](https://github.com/Conalog/sun_calculator)의 Go 버전입니다.

## 설치

[참고](https://devocean.sk.com/blog/techBoardDetail.do?ID=163883)

1. Github에 SSH Key를 등록하거나 Access Token을 발급받습니다
2. github private repo를 go env에 등록합니다
    ```bash
    $ go env -w GOPRIVATE=github.com/Conalog
    ```
3. git uri를 override합니다.
   module uri에 인증 정보가 포함되도록 github.com으로 들어왔을 때 instread option으로 인증 정보가 포함하도록 합니다.
   * Access token을 사용하는 경우
    ```bash
    $ git config --global url."https://${user}:${personal_access_token}@github.com/Conalog".insteadOf "https://github.com/Conalog"
   ```
   * SSH Key를 사용하는 경우
     1. ~/.ssh/config에 github을 추가합니다.
        ```
        Host github
            HostName github.com
            User git
            IdentityFile ~/.ssh/id_rsa
        ```
     2. git uri를 override합니다.
        ```bash
        $ git config --global url."github:Conalog".insteadOf "https://github.com/Conalog"
        ```
4. go mod init 된 디렉토리에서 `go get github.com/conalog/~`를 통해 모듈을 다운로드합니다.
